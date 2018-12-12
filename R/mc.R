#' Genarate metacells from TAD methylation
#' 
#' @param db cgdb object
#' @param tad_intervs intervals set of TADs
#' @param knn knn
#' @param k_expand k_expand
#' @param min_cluster_size min_cluster_size 
#' @param n_resamp n_resamp
#' @param seed seed 
#' @param K K
#' @param min_mc_size min_mc_size 
#' @param min_tad_cov min_tad_cov 
#' @param min_tad_cells min_tad_cells 
#' @param min_tads min_tads
#' @param tads_meth previously computed methylation of TADS (using get_tads_meth)
#' 
#' @return list with the following fields: mc_map, tads_clust, tads_meth, tads_meth_f, tads_mc  
#' 
#' @export
sc5mc.generate_metacells <- function(db, tad_intervs, knn=80, k_expand=10, min_cluster_size=20, n_resamp=1000, seed=NULL,  K=50, min_mc_size=30, min_tad_cov=10, min_tad_cells=50, min_tads=1000, tads_meth=NULL){
    if (is.null(tads_meth)){
        tads_meth <- get_tads_meth(db, tad_intervs=tad_intervs)    
    }
    
    tads_meth_f <- tads_meth %>% filter_tads(min_cov = min_tad_cov, min_cells=min_tad_cells, min_tads=min_tads)
    tads_clust <-  tads_meth_f %>% cluster_tads_knn(knn=knn, k_expand=k_expand, min_cluster_size=min_cluster_size, n_resamp=n_resamp, seed=seed, normalize=TRUE)
    tads_mc <- mc_from_coclust(tads_clust$coc, tads_clust$tads_m, K=K, min_mc_size=min_mc_size)

    mc_map <- tads_mc$node_clust %>% select(cell_id = node, mc = cluster) %>% mutate(mc = glue('mc{mc + 1}'))
    return(list(mc_map=mc_map, tads_clust=tads_clust, tads_meth=tads_meth, tads_meth_f=tads_meth_f, tads_mc=tads_mc))
}

#' @export
get_tads_meth <- function(db, tad_intervs){
    tads_meth <- db %>% summarise_intervals(tad_intervs)
    return(tads_meth)
}

#' @export
filter_tads <- function(tads_meth, min_cov, min_cells, min_tads=NULL){    
    tads_meth <- tads_meth %>% 
        group_by(chrom, start, end) %>% 
        filter(sum(cov >= min_cov) >= min_cells) 

    n_tads <- n_groups(tads_meth)    
    min_tads <- min_tads %||% (0.7 * n_tads)
    message(glue('min_tads: {min_tads}'))

    tads_meth <- tads_meth %>% 
        group_by(cell_id) %>% 
        filter(sum(cov >= min_cov) >= min_tads) %>%
        ungroup()

    return(tads_meth)
}

cluster_tads_knn <- function(tads_meth, knn, k_expand, min_cluster_size, n_resamp=500, seed=NULL, normalize=FALSE){
    if (normalize){
        tads_meth <- tads_meth %>% mutate(avg = meth / cov) %>% group_by(chrom, start, end) %>% mutate(avg = avg / sum(avg, na.rm=TRUE)) %>% ungroup()
    } else {
        tads_meth <- tads_meth %>% mutate(avg = meth / cov)   
    }
    tads_m <- tads_meth %>% select(chrom:end, cell_id, avg) %>% spread(cell_id, avg) %>% as.data.frame()
    tads_intervs <- tads_m %>% select(chrom, start, end)
    
    tads_m <- tads_m %>% 
        unite('tad', chrom:end) %>%
        column_to_rownames('tad') %>% 
        as.matrix()

    print(dim(tads_m))      

    cm <- tgs_cor(tads_m,  spearman=TRUE, pairwise.complete.obs=TRUE)
    cm_knn <- tgs_knn(cm, knn=knn*k_expand)    
    g <- tgs_graph(cm_knn, knn=knn, k_expand=k_expand)

    if (!is.null(seed)){
        set.seed(seed)
    }
    coc <- tgs_graph_cover_resample(graph=g, knn=knn, min_cluster_size=min_cluster_size, n_resamp=n_resamp)

    m_coc <- as.matrix(Matrix::sparseMatrix(as.numeric(coc$co_cluster$node1), as.numeric(coc$co_cluster$node2), x=coc$co_cluster$cnt, dims=c(length(coc$samples), length(coc$samples))))
    m_samp <- (coc$samples/mean(coc$samples)) %*% t(coc$samples)
    
    m_coc <- m_coc/m_samp
    rownames(m_coc) <- colnames(tads_m)
    colnames(m_coc) <- colnames(tads_m)
    outliers <- colnames(tads_m)[rowSums(m_coc)==0]   
    m_coc[lower.tri(m_coc)] <- t(m_coc)[lower.tri(m_coc)]    
    
    return(list(tads_m = tads_m, tads_intervs=tads_intervs, coclust=coc, m_coc = m_coc, outliers=outliers, cm=cm, cm_knn=cm_knn, g = g))
}

mc_from_coclust <- function(coc, mat, K, min_mc_size, alpha=2){
    edges <- coc$co_cluster
    filt_edges <- coclust_filt_by_k_deg(edges, K=K, alpha=alpha)
    message("filtered ", nrow(edges) - sum(filt_edges), " left with ", sum(filt_edges), " based on co-cluster imbalance")

    edges <- edges %>% 
        filter(filt_edges) %>%
        rlang::set_names(c("col1", "col2", "weight")) %>%
        mutate(weight = weight / max(weight)) %>% 
        filter(col1 != col2)

    edges <- bind_rows(edges, rename(edges, col1 = col2, col2 = col1))

    node_clust <- tgs_graph_cover(edges, min_mc_size)

    f_outlier <- (node_clust$cluster == 0)
    outliers <- colnames(mat)[node_clust$node[f_outlier]]

    mc <- as.integer(as.factor(node_clust$cluster))
    names(mc) <- colnames(mat)

    message("# of metacells: ", max(mc))

    return(list(mc = mc, outliers=outliers, node_clust=node_clust))
}

coclust_filt_by_k_deg <- function(coclust, K, alpha){
    edges <- coclust
    
    deg_wgt <- as.matrix(table(c(edges$node1, edges$node2), c(edges$cnt,edges$cnt)))
    deg_cum <- t(apply(deg_wgt, 1, function(x) cumsum(rev(x))))
    thresh_Kr <- rowSums(deg_cum > K)
    thresh_K <- rep(NA, length(levels(edges$node1)))
    names(thresh_K) <- levels(edges$node1)
    thresh_K[as.numeric(names(thresh_Kr))] <- thresh_Kr

    filt_edges <- thresh_K[edges$node1] < edges$cnt * alpha | 
                            thresh_K[edges$node2] < edges$cnt * alpha

    return(filt_edges)
}
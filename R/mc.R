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
#' @param rm_xy remove TADs from X and Y chromosome
#' @param normalize_tads normalize TAD methylation by dividing by the sum of each row
#' @param downsample downsample cells to have the same coverage
#' @param min_cov minimial coverage per cell (for downsampling). If NULL the value would be determined by the coverage of the cell in \code{ds_quantile} quantile
#' @param ds_quantile see \code{min_cov}. Quantile of the cell that determines the minimal coverage for downsampling.
#' @param norm_knn_mat matrix to use for normalizing orthogonal effects like cell cycle. Matrix should have a row per cell and a column per variable to normalize. For example, the matrix can have one column of coverage in early/coverage in late and another column of difference in methylation in early vs late regions.
#' @param norm_knn_k k parameter for normalization
#' 
#' 
#' @return list with the following fields: mc_map, tads_clust, tads_meth, tads_meth_f, tads_mc  
#' 
#' @export
sc5mc.generate_metacells <- function(db, tad_intervs, knn=80, k_expand=10, min_cluster_size=20, n_resamp=1000, seed=NULL,  K=50, min_mc_size=30, min_tad_cov=10, min_tad_cells=50, min_tads=1000, tads_meth=NULL, rm_xy = FALSE, normalize_tads = TRUE, downsample = FALSE, min_cov = NULL, ds_quantile = 0.1, norm_knn_mat = NULL, norm_knn_k = 20){

    if (is.null(tads_meth)){
        tads_meth <- get_tads_meth(db, tad_intervs=tad_intervs)    
    }

    if (rm_xy){
        tads_meth <- tads_meth %>% filter(!(chrom %in% c('chrX', 'chrY')))
    }
    
    tads_meth_f <- tads_meth %>% filter_tads(min_cov = min_tad_cov, min_cells=min_tad_cells, min_tads=min_tads)

    tads_mat <- tads_meth_f

    if (downsample){

        if (is.null(min_cov)){
            cells_n <- tads_meth_f %>% group_by(cell_id) %>% summarise(n = sum(cov))
            dsn <- quantile(cells_n$n, ds_quantile)
            message(glue('dsn = {dsn}'))            

            tads_meth_f <- tads_meth_f %>% filter(cell_id %in% cells_n$cell_id[cells_n$n >= dsn])
            tads_ds <- downsample_meth(tads_meth_f, dsn=dsn)    
        } else {            
            tads_meth_f <- tads_meth_f  %>% filter(cov >= min_cov)          
            tads_ds <- ds_loci_meth(tads_meth_f, min_cov)
            
            all_cell_num <- length(unique(tads_meth_f$cell_id) )
            tads_ds <- tads_ds %>% group_by(chrom, start, end) %>% filter(n() >= min_tad_cells)  %>% ungroup()
        }        
        
        tads_meth_f <- tads_ds
        tads_mat <- tads_ds
    }


    tads_clust <-  tads_mat %>% cluster_tads_knn(knn=knn, k_expand=k_expand, min_cluster_size=min_cluster_size, n_resamp=n_resamp, seed=seed, normalize=normalize_tads, norm_knn_mat = norm_knn_mat, norm_knn_k = norm_knn_k)
    tads_mc <- mc_from_coclust(tads_clust$coc, tads_clust$tads_m, K=K, min_mc_size=min_mc_size)

    mc_map <- tads_mc$node_clust %>% select(cell_id = node, mc = cluster) %>% mutate(mc = glue('mc{mc + 1}'))
    return(list(mc_map=mc_map, tads_clust=tads_clust, tads_meth=tads_meth, tads_meth_f=tads_meth_f, tads_mc=tads_mc, tads_ds = tads_ds, dsn = dsn))
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

#' @export
cluster_tads_knn <- function(tads_meth, knn, k_expand, min_cluster_size, n_resamp=500, seed=NULL, normalize=FALSE, norm_knn_mat = NULL, norm_knn_k = 20){    

    tads_m <- tads_meth %>% mutate(avg = meth / cov) %>% select(chrom:end, cell_id, avg) %>% spread(cell_id, avg) %>% as.data.frame()
    tads_intervs <- tads_m %>% select(chrom, start, end)
    
    tads_m <- tads_m %>% 
        unite('tad', chrom:end) %>%
        column_to_rownames('tad') %>% 
        as.matrix()

    print(dim(tads_m))      


    if (normalize){
        # tads_m <- tads_m - rowMeans(tads_m, na.rm=TRUE)         
        tads_m <- tads_m / rowSums(tads_m, na.rm=TRUE)    
    }    

    if (!is.null(norm_knn_mat)){
        norm_knn_mat[!is.finite(norm_knn_mat)] <- NA
        norm_knn_mat <- scale(norm_knn_mat, center = TRUE, scale = TRUE)                  
        
        cells <- intersect(rownames(norm_knn_mat), colnames(tads_m))
        norm_mat <- as.matrix(norm_knn_mat[cells, ])
        tads_m <- tads_m[, cells]

        cm <- tgs_cor(tads_m,  spearman=TRUE, pairwise.complete.obs=TRUE)
        
        knn_mat <- get_norm_knn_matrix(norm_mat, k = norm_knn_k)
        cm_norm <- normalize_knn(cm, knn_mat)                 
        cm_knn <- tgs_knn(cm_norm, knn=knn*k_expand)      
    } else {               
        cm <- tgs_cor(tads_m,  spearman=TRUE, pairwise.complete.obs=TRUE)
        cm_norm <- NULL 
        norm_mat <- NULL
        cm_knn <- tgs_knn(cm, knn=knn*k_expand)    
    }    

    
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
    
    return(list(tads_m = tads_m, tads_intervs=tads_intervs, coclust=coc, m_coc = m_coc, outliers=outliers, cm=cm, cm_knn=cm_knn, g = g, cm_norm = cm_norm, norm_knn_mat = norm_mat, norm_knn_k = norm_knn_k))
}

#' @export
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

#' @export
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

ds_loci_meth <- function(df, min_cov){
    df <- df %>% filter(cov >= min_cov)
    df %>% mutate(cov = min_cov, meth = map2_int(df$cov, df$meth, ~ sum(sample(1:.x, min_cov) <= .y)) )
} 

normalize_knn <- function(raw, knn_mat) {
    raw <- raw[colnames(knn_mat), colnames(knn_mat)]
    raw_filt_na <- raw
    raw_filt_na[is.na(raw)] <- 0
    met_exp <- as.matrix(raw_filt_na) %*% as.matrix(t(knn_mat))
    # met_exp[lower.tri(met_exp)] <- t(met_exp)[lower.tri(met_exp)]
    not_na_mat <- !is.na(raw)
    met_exp_n <- as.matrix(not_na_mat) %*% as.matrix(t(knn_mat))
    # met_exp_n[lower.tri(met_exp_n)] <- t(met_exp_n)[lower.tri(met_exp_n)]
    met_exp_norm <- met_exp / met_exp_n
    met_oe <- raw - met_exp_norm
    met_oe[lower.tri(met_oe)] <- t(met_oe)[lower.tri(met_oe)]   
    return(met_oe)
}

get_norm_knn_matrix <- function(feats, k) {
    dist_mat <- tgs_dist(feats)
    knn_df <- tgs_knn(100 - as.matrix(dist_mat), k)
    knn_df <- knn_df %>% mutate(col1 = factor(col1), col2 = factor(col2, levels = levels(col1)))

    knn_mat <- sparseMatrix(as.numeric(knn_df$col1), as.numeric(knn_df$col2), x = 1)
    rownames(knn_mat) <- levels(knn_df$col1)
    colnames(knn_mat) <- levels(knn_df$col2)
    return(knn_mat)
}

#' Cluster cells
#' 
#' @param clustering_distance distance measure used for clustering the correlations
#' @param clustering_method  clustering method used for clustering the correlations. Accepts the same values as ‘hclust’.
#' @param plot plot the clustered correlation matrix
#' @param ... additional parameters for sc5mc.plot_cor_mat
#' @inheritParams sc5mc.calc_pdiff
#' @inheritParams sc5mc.calc_pdiff_cor
#' @inheritParams sc5mc.plot_cor_mat
#' 
#' @return list with the following components: 
#' \describe{
#'   \item{cell_comp:}{pairwise comparison of cell CpGs.}
#'   \item{cell_cor:}{pairwise cell correlations.}
#'   \item{cell_cor_mat:}{cells correlation matrix.}
#'   \item{cell_cor_hclust_rows:}{hclust object used to cluster the rows.}
#'   \item{cell_cor_hclust_cols:}{hclust object used to cluster the columns.}
#' }
#' If plot is TRUE the smat object is returned invisibly
#' 
#' @export
sc5mc.cluster_cells <- function(smat, min_cells=10, min_cgs=100, pairwise.complete.obs=TRUE, spearman=TRUE, clustering_distance = 'euclidean',  clustering_method = 'ward.D2', intervs=NULL, cols=NULL, samp_data=NULL, plot=FALSE, ...){

    smat <- smat.cluster_cells(smat=smat, min_cells=min_cells, min_cgs=min_cgs, pairwise.complete.obs=pairwise.complete.obs, spearman=spearman, clustering_distance=clustering_distance, clustering_method=clustering_method, intervs=intervs, cols=cols, samp_data=samp_data, plot=plot, ...)
  
    res <- list(cell_comp = smat$cell_comp,
                cell_cor = smat$cell_cor, 
                cell_cor_mat = smat$cell_cor_mat, 
                cell_cor_hclust_rows = smat$cell_cor_hclust_rows, 
                cell_cor_hclust_cols = smat$cell_cor_hclust_cols)

    return(res)
}

#' Cluster cells
#' 
#' @inheritParams sc5mc.cluster_cells
#' 
#' @return smat object with the following fields: 
#' \describe{
#'   \item{cell_comp:}{pairwise comparison of cell CpGs.}
#'   \item{cell_cor:}{pairwise cell correlations.}
#'   \item{cell_cor_mat:}{cells correlation matrix.}
#'   \item{cell_cor_hclust_rows:}{hclust object used to cluster the rows.}
#'   \item{cell_cor_hclust_cols:}{hclust object used to cluster the columns.}
#' }
#' If plot is TRUE the smat object is returned invisibly
#' 
#' @export
smat.cluster_cells <- function(smat, min_cells=10, min_cgs=100, pairwise.complete.obs=TRUE, spearman=TRUE, clustering_method = 'ward.D2', intervs=NULL, cols=NULL, samp_data=NULL, plot=FALSE, ...){
    smat <- smat.calc_pdiff_cor(smat=smat, min_cells=min_cells, min_cgs=min_cgs, pairwise.complete.obs=pairwise.complete.obs, spearman=spearman, intervs=intervs, cols=cols, samp_data=samp_data)

    smat$cell_cor_mat <- smat[['cell_cor']] %>% filter(cell1 != cell2) %>% spread(cell2, corr) %>% select(-cell1) %>% as.matrix()
    mat <- smat$cell_cor_mat

    message('clustering...')
    smat$cell_cor_hclust_cols <- tgstat::tgs_dist(t(mat)) %>% hclust(method = clustering_method) %>% dendsort::dendsort()
    smat$cell_cor_hclust_rows <- tgstat::tgs_dist(mat) %>% hclust(method = clustering_method) %>% dendsort::dendsort()

    if (plot){
        sc5mc.plot_cor_mat(smat$cell_cor, cluster_rows=smat$cell_cor_hclust_rows, cluster_cols=smat$cell_cor_hclust_cols, ...)
        invisible(smat)
    }
    return(smat)
}

#' calculate pairwise statistics of cells
#'
#' @param smat smat object
#'
#' @param min_cgs minimal number of cpgs covered in both cells
#' @param intervs intervals set to filter smat by before calculating. if NULL
#' all CpGs would be used.
#' @param cols vector of names of cells to use. If NULL all cells would be used.
#' @param samp_data data frame with additional annotation of the cell.
#' must have 'track' field with the cell name.
#'
#' @export
sc5mc.calc_pdiff <- function(smat, min_cgs=100, intervs=NULL, cols=NULL, samp_data=NULL){
    if (!is.null(intervs) || !is.null(cols)){
        smat <- smat.filter(smat, intervs=intervs, cols=cols)
    }

    calc_mat_stat <- function(smat, mat1_name, mat2_name, column_name){
        if (mat1_name == mat2_name){
            stat <- crossprod(smat[[mat1_name]])
        } else {
            stat <- crossprod(smat[[mat1_name]], smat[[mat2_name]])
        }
        stat <- stat %>% as.matrix() %>% reshape2::melt() %>% purrr::set_names(c('cell1', 'cell2', column_name))
        return(stat)
    }

    stats_tab <- tribble(
       ~mat1,    ~mat2,        ~col,
       'meth',   'meth',       'n11',
       'unmeth', 'unmeth',     'n00',
       'meth',   'unmeth',     'n10',
       'unmeth', 'meth',       'n01',
       'cov',    'cov',        'ntot' )

    cell_comp <- plyr::alply(stats_tab, 1, function(x) calc_mat_stat(smat, x$mat1, x$mat2, x$col), .parallel = TRUE)

    cell_comp <- cell_comp %>%  reduce(left_join, by=c('cell1', 'cell2')) %>% tbl_df

    cell_comp <- cell_comp %>%
        mutate(n = n11 + n00 + n10 + n01, psame = (n00 + n11) / n, pdiff=1 - psame) %>%
        filter(cell1 != cell2, n >= min_cgs)

    if (!is.null(samp_data)){
        cell_comp <- cell_comp %>% left_join(samp_data %>% select(cell1 = track, lib), by='cell1') %>% mutate(cell1 = lib) %>% select(-lib)
        cell_comp <- cell_comp  %>% left_join(samp_data %>% select(cell2 = track, lib), by='cell2') %>% mutate(cell2 = lib) %>% select(-lib)
    }

    return(cell_comp)
}

#' calculate pairwise statistics of cells
#'
#' @param smat smat object
#'
#' @inheritParams sc5mc.calc_pdiff
#' 
#' @return smat object with cell_comp field 
#'
#' @export
smat.calc_pdiff <- function(smat, min_cgs=100, intervs=NULL, cols=NULL, samp_data=NULL){
    message('calculating pairwise statistics of cells...')
    smat$cell_comp <- sc5mc.calc_pdiff(smat, min_cgs=min_cgs, intervs=intervs, cols=cols, samp_data=samp_data)
    return(smat)
}

#' Calculate pairwise correlation of cells pdiff
#' 
#' @param cell_comp Output of \code{sc5mc.calc_pdiff}: data frame with pairs of cells (cell1, cell2) and 'pdiff' column with the fraction of different CpGs.
#'
#' @param min_cells minimal number of cells with pdiff statistic in both cells
#' @param pairwise.complete.obs if TRUE: similar to use = 'pairwise.complete.obs' in \code{cor}. if FALSE: similar to use = 'everything'
#' @param spearman if 'TRUE' Spearman correlation is computed, otherwise Pearson
#' @param cells vector with cell names of cells to operate on. If NULL calculation would be 
#' done on all cells.
#' 
#'
#' @export
sc5mc.calc_pdiff_cor <- function(cell_comp, min_cells, pairwise.complete.obs=TRUE, spearman=TRUE, cells=NULL){
    cell_comp <- cell_comp %>% select(cell1, cell2, pdiff)

    if (!is.null(cells)){
        cell_comp <- cell_comp %>% filter(cell1 %in% cells, cell2 %in% cells)
    }

    good_cells <- cell_comp %>%
        group_by(cell1) %>%
        mutate(val_num = n()) %>%
        filter(val_num >= min_cells) %>%
        .$cell1 %>%
        unique()

    cell_comp <- cell_comp %>%
        filter(cell1 %in% good_cells, cell2 %in% good_cells)

    cell_comp <- cell_comp  %>%
        spread(cell2, pdiff) %>%
        .[,-1] %>%
        as.matrix()

    cell_cor <- tgstat::tgs_cor(cell_comp, pairwise.complete.obs = pairwise.complete.obs, spearman = spearman) %>%
        reshape2::melt() %>%
        rename(cell1=Var1, cell2=Var2, corr=value) %>%
        tbl_df
    return(cell_cor)
}

#' calculate pairwise correlation of cells pdiff
#'
#' @param smat smat object
#'
#' @inheritParams smat.calc_pdiff
#' @inheritParams sc5mc.calc_pdiff_cor
#' 
#' @return smat object with cell_cor field 
#'
#' @export
smat.calc_pdiff_cor <- function(smat, min_cells=10, min_cgs=100, pairwise.complete.obs=TRUE, spearman=TRUE, intervs=NULL, cols=NULL, samp_data=NULL){
    if (!('cell_comp' %in% names(smat))){        
        smat  <- smat.calc_pdiff(smat, min_cgs=min_cgs, intervs=intervs, cols=cols, samp_data=samp_data)
    }    

    message('calculating pairwise correlations...')    
    smat$cell_cor <- sc5mc.calc_pdiff_cor(smat$cell_comp, min_cells=min_cells, pairwise.complete.obs=pairwise.complete.obs, spearman=spearman)
    return(smat)
}


#' calculate pairwise statistics of CpGs
#' 
#' @param smat smat object
#'
#' @param min_cells minimal number of cells
#' @param intervs 
#' @param cols
#'
#' @export
sc5mc.calc_pdiff_rows <- function(smat, min_cells=5, intervs=NULL, cols=NULL){
    if (!is.null(intervs) || !is.null(cols)){
        smat <- smat.filter(smat, intervs=intervs, cols=cols)
    }

    calc_mat_stat <- function(smat, mat1_name, mat2_name, column_name){
        if (mat1_name == mat2_name){
            stat <- tcrossprod(smat[[mat1_name]])
        } else {
            stat <- tcrossprod(smat[[mat1_name]], smat[[mat2_name]])
        }
        stat <- stat %>% as.matrix() %>% reshape2::melt() %>% purrr::set_names(c('interval1', 'interval2', column_name))
        return(stat)
    }

    stats_tab <- tribble(
        ~mat1, ~mat2, ~col,
        'meth',   'meth',       'n11',
        'unmeth', 'unmeth',     'n00',
        'meth',   'unmeth',     'n10',
        'unmeth', 'meth',       'n01',
        'cov',    'cov',        'ntot'
    )

    interval_comp <- plyr::alply(stats_tab, 1, function(x) calc_mat_stat(smat, x$mat1, x$mat2, x$col), .parallel=T)


    interval_comp <- interval_comp %>%  reduce(left_join, by=c('interval1', 'interval2')) %>% tbl_df

    interval_comp <- interval_comp %>%
        mutate(n = n11 + n00 + n10 + n01, psame = (n00 + n11) / n, pdiff=1 - psame) %>%
        filter(interval1 != interval2, n >= min_cells)


    return(interval_comp)
}

#' calculate pairwise correlation of CpGs pdiff
#' 
#' @param row_comp
#'
#' @param min_intervals min_intervals
#' @param spearman spearman
#' @param pairwise.complete.obs pairwise.complete.obs
#' @param intervals intervals
#'
#' @export
sc5mcs.calc_pdiff_cor_rows <- function(row_comp, min_intervals, pairwise.complete.obs=TRUE, spearman=TRUE, intervals=NULL){
    row_comp <- row_comp %>% rename(cell1 = interval1, cell2 = interval2)
    row_cor <- sc5mc.calc_pdiff_cor(row_comp, min_cells=min_intervals, pairwise.complete.obs=pairwise.complete.obs, spearman=spearman, cells=intervals)
    row_cor <- row_cor %>% rename(interval1 = cell1, interval2 = cell2)
    return(row_cor)
}


#' plot cell cell correlation matrix
#' 
#' @param cell_cor cell cell correlation (output of sc5mc.calc_pdiff_cor)
#'
#' @param row_ord specific order for the rows
#' @param col_ord specific order for the columns
#' @param show_colnames show column names
#' @param show_rownames show row names
#' @param color_pal color pallete
#' @param ... other parameters of pheatmap1
#' @inheritParams gpatterns::pheatmap1
#'
#' @export
sc5mc.plot_cor_mat <- function(cell_cor,
                               row_ord = NULL,
                               col_ord = NULL,
                               show_colnames = FALSE,
                               show_rownames = FALSE,
                               breaks = NULL,
                               color_pal = NULL,
                               ...) {
    if (is.null(color_pal)){
        color_pal <- c('blue', 'white', '#ffffd4','#fed98e','#fe9929','#d95f0e','#993404', 'black')
        color_breaks <-  quantile(cell_cor[cell_cor$cell1 != cell_cor$cell2, ]$corr, seq(0,1,1/(length(color_pal) - 1)), na.rm=T)
        pallete <- build_pallette(data.frame(point =color_breaks, color =color_pal), 1000)
    }

    if (is.null(breaks)){
        color_breaks <-  quantile(cell_cor[cell_cor$cell1 != cell_cor$cell2, ]$corr, seq(0,1,1/(length(color_pal) - 1)), na.rm=T)
        zlim <- c(min(color_breaks), max(color_breaks))
        breaks <- seq(zlim[1], zlim[2], length.out=1001)
    }

    cell_cor <- cell_cor %>% filter(cell1 != cell2) %>% spread(cell2, corr)
    if (!is.null(row_ord)){
        cell_cor <- cell_cor[row_ord, ]
    }

    if (!is.null(col_ord)){
        cell_cor <- cell_cor[, c(1, col_ord + 1)]
    }

    res <- cell_cor %>% pheatmap1(color=pallete, method='ward.D2', breaks=breaks, show_colnames=show_colnames, show_rownames=show_rownames, ...)
    invisible(res)
}

#' plot cell cell correlation matrix
#' 
#' @param smat smat object
#' @param title show title
#' @inheritParams sc5mc.plot_cor_mat
#'
#' @export
smat.plot_cor_mat <- function(smat,
                               row_ord = NULL,
                               col_ord = NULL,
                               show_colnames = FALSE,
                               show_rownames = FALSE,
                               breaks = NULL,
                               color_pal = NULL,    
                               title = TRUE,                           
                               ...) {
    if (!('cell_cor' %in% names(smat))){
        stop('No "cell_cor" field. Please run smat.cluster_cells() first')
    }

    if ('cell_cor_hclust_rows' %in% names(smat)){
        cluster_rows <- smat$cell_cor_hclust_rows
    } else {
        cluster_rows <- FALSE
    }

    if ('cell_cor_hclust_cols' %in% names(smat)){
        cluster_cols <- smat$cell_cor_hclust_cols
    } else {
        cluster_cols <- FALSE
    }

    if (title){
        if (has_name(smat, 'name')){
            name_str <- glue('{smat$name}\n')
        } else {
            name_str <- ''
        }
        cell_num <- length(unique(smat_f$cell_cor$cell2))
        main <- glue('{name_str}{comify(cell_num)} cells')
    } else {
        main <- NA
    }

    sc5mc.plot_cor_mat(smat$cell_cor, row_ord=row_ord, col_ord=col_ord, show_rownames=show_colnames, breaks=breaks, color_pal=color_pal, cluster_rows=cluster_rows, cluster_cols=cluster_cols, main=main,  ...)
}


shuffle_mat <- function(m_meth, m_cov, m_avg, n_shuff=1e3){
    m_unmeth <- m_cov - m_meth

    m_shuff <- shuffle_mat_marginals(m_meth, which(m_meth > 0, arr.ind=TRUE), n_shuff )

    m_cov_shuff <- (m_shuff + m_unmeth)
    m_avg_shuff <- m_shuff / (m_shuff + m_unmeth)

    suppressWarnings(cr <- cor.test(as_vector(m_avg_shuff[m_cov > 0]), as_vector(m_avg[m_cov > 0]), method='spearman', na.rm=T)$estimate            )
    stopifnot(all(colSums(m_shuff) == colSums(m_meth)))
    stopifnot(all(colSums(m_shuff) == colSums(m_meth)))
    message(glue('cor = {cr}'))
    return(list(cov = m_cov_shuff, meth = m_shuff, avg = m_avg_shuff))
}



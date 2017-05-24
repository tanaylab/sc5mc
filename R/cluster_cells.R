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
        stat <- stat %>% as.matrix() %>% reshape2::melt() %>% set_names(c('cell1', 'cell2', column_name))
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

    cell_comp <- plyr::alply(stats_tab, 1, function(x) calc_mat_stat(smat, x$mat1, x$mat2, x$col), .parallel=T)

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
        stat <- stat %>% as.matrix() %>% reshape2::melt() %>% set_names(c('interval1', 'interval2', column_name))
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

#' @export
sc5mcs.calc_pdiff_cor_rows <- function(row_comp, min_intervals, use='pairwise.complete.obs', method='spearman', intervals=NULL){
    row_comp <- row_comp %>% rename(cell1 = interval1, cell2 = interval2)
    row_cor <- sc5mc.calc_pdiff_cor(row_comp, min_cells=min_intervals, use=use, method=method, cells=intervals)
    row_cor <- row_cor %>% rename(interval1 = cell1, interval2 = cell2)
    return(row_cor)
}

#' @export
sc5mc.calc_pdiff_cor <- function(cell_comp, min_cells, use='pairwise.complete.obs', method='spearman', cells=NULL){
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

    cell_cor <- cor(cell_comp, use=use, method=method) %>% 
        reshape2::melt() %>% 
        rename(cell1=Var1, cell2=Var2, corr=value) %>% 
        tbl_df
    return(cell_cor)
}

#' @export
sc5mc.plot_cor_mat <- function(cell_cor, min_vals_row=100, min_vals_col=100, row_ord=NULL, col_ord=NULL,  show_colnames=FALSE, show_rownames=FALSE, breaks=NULL, color_pal=NULL, ...){

    if (is.null(color_pal)){
        color_pal <- c('blue', 'white', '#ffffd4','#fed98e','#fe9929','#d95f0e','#993404', 'black')   
        breaks1 <-  quantile(cell_cor[cell_cor$cell1 != cell_cor$cell2, ]$corr, seq(0,1,1/(length(color_pal) - 1)), na.rm=T)
        pallete <- build_pallette(data.frame(point =breaks1, color =color_pal), 1000)           
    }

    if (is.null(breaks)){
        breaks1 <-  quantile(cell_cor[cell_cor$cell1 != cell_cor$cell2, ]$corr, seq(0,1,1/(length(color_pal) - 1)), na.rm=T)
        zlim <- c(min(breaks1), max(breaks1))
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



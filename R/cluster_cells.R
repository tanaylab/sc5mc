#' @export
calc_pdiff <- function(smat, samp_data=NULL, min_cgs=100){
    n11 <- crossprod(smat$meth) %>% as.matrix %>% melt %>% rename(cell1=Var1, cell2=Var2, n11=value)
    n00 <- crossprod(smat$unmeth) %>% as.matrix %>% melt %>% rename(cell1=Var1, cell2=Var2, n00=value)
    n10 <- crossprod(smat$meth, smat$unmeth) %>% as.matrix %>% melt %>% rename(cell1=Var1, cell2=Var2, n10=value)
    n01 <- crossprod(smat$unmeth, smat$meth) %>% as.matrix %>% melt %>% rename(cell1=Var1, cell2=Var2, n01=value)
    ntot <- crossprod(smat$cov) %>% as.matrix %>% melt %>% rename(cell1=Var1, cell2=Var2, ntot=value)

    cell_comp <- list(n11, n00, n10, n01, ntot) %>% reduce(left_join) %>% tbl_df  
    cell_comp <- cell_comp %>% mutate(n = n11 + n00 + n10 + n01, psame = (n00 + n11) / n, pdiff=1 - psame) %>% filter(cell1 != cell2, n >= min_cgs) #%>% select(cell1, cell2, psame, pdiff)    
    
    if (!is.null(samp_data)){
        cell_comp <- cell_comp %>% left_join(samp_data %>% rename(cell1 = track)) %>% mutate(cell1 = samp) %>% select(-group, -samp) %>% left_join(samp_data %>% rename(cell2 = track)) %>% mutate(cell2 = samp) %>% select(-group, -samp)        
    }
    
    return(cell_comp)
}

#' @export
calc_pdiff_cor <- function(cell_comp, min_cells){
    good_cells <- cell_comp %>% group_by(cell1) %>% mutate(val_num = n()) %>% filter(val_num >= min_cells) %>% .$cell1 %>% unique
    cell_comp <- cell_comp %>% filter(cell1 %in% good_cells, cell2 %in% good_cells)
    cells_cor <- cell_comp  %>% select(-psame) %>% spread(cell2, pdiff) %>% .[,-1] %>% as.matrix %>% cor(use='pairwise.complete.obs',  method='spearman') %>% melt %>% rename(cell1=Var1, cell2=Var2, corr=value) %>% tbl_df
    return(cells_cor)
}

#' @export
plot_cor_mat <- function(mat, min_vals_row=100, min_vals_col=100, row_ord=NULL, col_ord=NULL, show_rownames=F, show_colnames=F, ...){
    annots <- samp_data %>% separate(group, c('plate', 'type')) %>% select(samp, plate, type)
    annot_colors <- data.frame(type = 'type', variable=c('p0', 'p1'), color=c('green', 'red')) %>% bind_rows(data.frame(type = 'plate', variable=c('plate1', 'plate7'), color=c('purple', 'yellow')))

    shades <- c('blue', 'white', '#ffffd4','#fed98e','#fe9929','#d95f0e','#993404', 'black')
    breaks <-  quantile(mat[mat$cell1 != mat$cell2, ]$corr, seq(0,1,1/(length(shades) - 1)), na.rm=T)
    pallete <- build_pallette(data.frame(point =breaks, color =shades), 1000)   
    zlim <- c(min(breaks), max(breaks))
    breaks <- seq(zlim[1], zlim[2], length.out=1001)
    
    # d_rows <- mat %>% filter(cell1 != cell2) %>% mutate(corr = log2(1 + 1 - corr)) %>% spread(cell2, corr) %>% select(-cell1) %>% as.dist
    # d_cols <- mat %>% filter(cell1 != cell2) %>% mutate(corr = log2(1 + 1 - corr)) %>% spread(cell2, corr) %>% select(-cell1) %>% t %>% as.dist   
    mat <- mat %>% filter(cell1 != cell2) %>% spread(cell2, corr) 
    if (!is.null(row_ord)){
        mat <- mat[row_ord, ]
    }

    if (!is.null(col_ord)){
        mat <- mat[, c(1, col_ord + 1)]
    }

    res <- mat %>% pheatmap1(show_rownames=show_rownames, show_colnames=show_colnames, annotation=annots, annotation_col=c('type', 'plate'), annotation_colors=annot_colors, color=pallete, method='ward.D2', breaks=breaks, ...)   
    invisible(res)
}



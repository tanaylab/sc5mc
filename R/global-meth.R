
#' Calculate average methylation per cell given CpG content
#' 
#' @param db cgdb object
#' @param min_cov minimal number of CpGs for average methylation calculation
#' @param breaks breaks of CpG content
#' @param downsample downsample the cells. Note that cells that do not have enough covered CpGs would be removed. 
#' @param downsample_n number of CpGs per each cell. If NULL - the number of CpGs would be determined based on the 0.1 quantile of cell coverage in each strata of CpG content. 
#' 
#' @return cell metadata of \code{db} with additional columns containing the average methylation in 
#' different CpG content regimes. If \code{downsample} is TRUE an additional column 'downsample_n' with the number of cpgs per strata of CpG content. 
#' 
#' @export
get_cgc_trend <- function(db, min_cov=300, breaks=c(seq(0,0.08,0.01), 0.15, 1), downsample=FALSE, downsample_n=NULL){
    cgc <- db %>% mutate_cpgs(cg_cont = cut(cg500, breaks=breaks)) %>% group_by_cpgs(cg_cont, add=TRUE) %>% summarise() %>% filter(cov >= min_cov)

    if (downsample){
        if (is.null(downsample_n)){
            dsn <- cgc %>% group_by(cg_cont) %>% summarise(downsample_n = quantile(cov, 0.1))
        } else {
            dsn <- cgc %>% group_by(cg_cont) %>% summarise(downsample_n = downsample_n)
        }
        cgc <- cgc %>% left_join(dsn) %>% plyr::ddply(plyr::.(cg_cont), function(x) downsample_cpgs(select(x, cell_id, cov, meth), n=x$downsample_n), .parallel=TRUE) %>% as_tibble()         
    }

    d <- cgc %>% mutate(avg = meth / cov)  %>% left_join(db@cells) %>% select(-cov, -meth) %>% spread(cg_cont, avg) 

    if (downsample){
        dsn <- dsn %>% mutate(cg_cont = paste0(cg_cont, '_dsn')) %>% spread(cg_cont, downsample_n)
        d <- cbind(d, dsn) %>% as_tibble()
    }    

    return(d)
}


plot_cg_cont_trend <- function(d, x='`(0.01,0.02]`', y='`(0.02,0.03]`', xlab=NULL, ylab=NULL, color_column=NULL, point_size=0.5, plot_sd=TRUE, sd_per_group=FALSE){       
    if (length(groups(d)) > 0){     
        color_column <- color_column %||% 'type'
        d <- d %>% unite_(color_column, groups(d))
    }

    breaks_to_label <- function(str){
        stringr::str_split(str, ',')[[1]] %>% gsub('`', '', .) %>% gsub('\\(', '', .) %>% gsub('\\[', '', .)  %>% gsub('\\)', '', .) %>% gsub('\\]', '', .) %>% gsub(' ', '', .) %>% as.numeric() %>% scales::percent(accuracy=1) %>% paste(collapse = '-')
    }

    xlab <- xlab %||% breaks_to_label(x)
    ylab <- ylab %||% breaks_to_label(y)

    x_col <- gsub("`", "", x)
    x_dsn_col <- glue('{x_col}_dsn')
    if (x_dsn_col %in% colnames(d)){
        x_dsn <- mean(d[[x_dsn_col]], na.rm=TRUE)
        xlab <- glue('{xlab} ({scales::comma(x_dsn)})')
    }

    y_col <- gsub("`", "", y)
    y_dsn_col <- glue('{y_col}_dsn')
    if (y_dsn_col %in% colnames(d)){
        y_dsn <- mean(d[[y_dsn_col]], na.rm=TRUE)
        ylab <- glue('{ylab} ({scales::comma(y_dsn)})')
    }
    
    p <- d %>% ggplot(aes_string(x=x, y=y, color=color_column)) + geom_point(size=0.5, shape=19, alpha=0.5) + xlab(xlab) + ylab(ylab) + ggsci::scale_color_lancet()

    if (plot_sd && x_dsn_col %in% colnames(d)){
        if (sd_per_group){            
            d_sd_lines <- d %>% select(cell_id, !! x_col, !! x_dsn_col, !! color_column) %>% group_by(!!! rlang::syms(color_column)) %>% do(get_meth_expected_sd(.[[x_col]], dsn=x_dsn)) %>% select(!! color_column, sd_upper, sd_lower) %>% gather('x', 'xintercept', - !! color_column)
            p <- p + geom_vline(data = d_sd_lines, aes_string(color=color_column, xintercept='xintercept'), linetype='dashed', alpha=0.5)    
        } else {
            exp_sd <- get_meth_expected_sd(d[[x_col]], dsn=x_dsn)
            p <- p + geom_vline(xintercept = c(exp_sd$sd_upper, exp_sd$sd_lower), linetype='dashed', color='black')
        }
    }

    if (plot_sd && y_dsn_col %in% colnames(d)){
        if (sd_per_group){
            d_sd_lines <- d %>% select(cell_id, !! y_col, !! y_dsn_col, !! color_column) %>% group_by(!!! rlang::syms(color_column)) %>% do(get_meth_expected_sd(.[[y_col]], dsn=y_dsn)) %>% select(!! color_column, sd_upper, sd_lower) %>% gather('x', 'yintercept', - !! color_column)
            p <- p + geom_hline(data = d_sd_lines, aes_string(color=color_column, yintercept='yintercept'), linetype='dashed', alpha=0.5)
        } else {
            exp_sd <- get_meth_expected_sd(d[[y_col]], dsn=y_dsn)
            p <- p + geom_hline(yintercept = c(exp_sd$sd_upper, exp_sd$sd_lower), linetype='dashed', color='black')            
        }        
    }

    return(p)
}

get_meth_expected_sd <- function(meth, dsn){
    p <- mean(meth, na.rm=TRUE)
    exp_sd <- sqrt(p*(1-p)*dsn)
    return(tibble(p = p, sd = exp_sd, sd_upper = p + exp_sd / dsn, sd_lower = p - exp_sd / dsn))
}

plot_cg_cont_trend_groups <- function(db, min_cov=300, breaks=seq(0,0.15,0.005), xlab='CpG content', ylab='Average methylation', ylim=c(0,1)){
    group_column <- cell_groups(db)

    d <- db %>% 
        mutate_cpgs(cg_cont = cut(cg500, breaks=breaks, include.lowest=TRUE)) %>%
        group_by_cpgs(cg_cont) %>% 
        summarise() %>%
        filter(cov >= min_cov) %>%
        filter(!is.na(cg_cont)) %>%
        mutate(avg = meth / cov)    

    p <- d %>% ggplot(aes_string(x='cg_cont', y='avg', color=group_column, group=group_column)) + 
            geom_line() + 
            ggsci::scale_color_lancet() + 
            xlab(xlab) + 
            ylab(ylab) + 
            ylim(ylim)

    return(p)
}


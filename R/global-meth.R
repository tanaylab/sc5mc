
#' Calculate average methylation per cell given CpG content
#' 
#' @param db cgdb object
#' @param min_cov minimal number of CpGs for average methylation calculation
#' @param breaks breaks of CpG content#' 
#' 
#' @return cell metadata of \code{db} with additional columns containing the average methylation in 
#' different CpG content regimes
#' 
#' @export
get_cgc_trend <- function(db, min_cov=300, breaks=c(seq(0,0.08,0.01), 0.15, 1)){
    cgc <- db %>% mutate_cpgs(cg_cont = cut(cg500, breaks=breaks)) %>% group_by_cpgs(cg_cont, add=TRUE) %>% summarise()

    d <- cgc %>% filter(cov >= min_cov) %>% mutate(avg = meth / cov)  %>% left_join(db@cells)  %>% select(-cov, -meth) %>% spread(cg_cont, avg) 
    return(d)
}

plot_cg_cont_trend <- function(d, x='`(0.01,0.02]`', y='`(0.02,0.03]`', xlab=NULL, ylab=NULL, color_column=NULL){		
	if (length(groups(d)) > 0){		
		color_column <- color_column %||% 'type'
		d <- d %>% unite_(color_column, groups(d))
	}

	breaks_to_label <- function(str){
		stringr::str_split(str, ',')[[1]] %>% gsub('`', '', .) %>% gsub('\\(', '', .) %>% gsub('\\[', '', .)  %>% gsub('\\)', '', .) %>% gsub('\\]', '', .) %>% gsub(' ', '', .) %>% as.numeric() %>% scales::percent() %>% paste(collapse = '-')
	}

	xlab <- xlab %||% breaks_to_label(x)
	ylab <- ylab %||% breaks_to_label(y)
	
	p <- d %>% ggplot(aes_string(x=x, y=y, color=color_column)) + geom_point(size=0.5, shape=19, alpha=0.5) + xlab(xlab) + ylab(ylab) + ggsci::scale_color_lancet()
	return(p)
}


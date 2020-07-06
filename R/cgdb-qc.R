#' Plot QC statistics of cgdb object
#' 
#' @param db cgdb object
#' @param ofn output file
#' @param min_cgc_cov minimal coverage for average methylation per CpG content estimation
#' @param return_all_figs return a list of all the figures
#' @param downsample downsample cells when plotting global trends
#' @param ... additional parameters for cowplot::save_plot
#' @inheritParams cowplot::save_plot
#' 
#' @return if \code{return_all_plots} - returns a list of all the sub-figures. otherwise returns a large plot with all the figures combined. 
#' 
#' 
#' @export
qc_plot <- function(db, ofn=NULL, base_height=12, min_cgc_cov=300, return_all_figs = FALSE, downsample=TRUE, ...){
	stats <- db@cells

	plots <- list()

	plots$p_seq <- sc5mc.plot_reads_per_cpg(db, log_scale=TRUE) + guides(color=FALSE)
	
	loginfo('calculating average methylation per CpG content')
	cgc_meth <- db %>% ungroup_cells() %>% ungroup_cpgs() %>% get_cgc_trend(min_cov=min_cgc_cov, breaks=c(seq(0,0.08,0.01), 0.12, 1), downsample=downsample)
	if (length(groups(stats)) > 0){		
		cgc_meth <- cgc_meth %>% group_by(!!! groups(stats))	
	}	

	loginfo('calculating CpG marginals')
	plots$p_cpg_cov <- sc5mc.plot_cpg_marginals_bars(db)

	loginfo('calculating cell marginals')
	plots$p_cell_cov <- sc5mc.plot_cell_marginals_bars(db)

	plots$p_cgc_low <- plot_cg_cont_trend(cgc_meth, x='`(0.01,0.02]`', y='`(0.02,0.03]`') + coord_cartesian(xlim=c(0,1), ylim=c(0,1))
	plots$p_cgc_low_high <- plot_cg_cont_trend(cgc_meth, x='`(0.01,0.02]`', y='`(0.08,0.12]`') + guides(color=FALSE) + coord_cartesian(xlim=c(0,1), ylim=c(0,0.3))

	high_cgc_meth <- db %>% ungroup_cells() %>% ungroup_cpgs() %>% get_cgc_trend(min_cov=min_cgc_cov, breaks=c(seq(0.08,0.1,0.01), 1), downsample=downsample) 
	if (length(groups(stats)) > 0){		
		high_cgc_meth <- high_cgc_meth %>% group_by(!!! groups(stats))	
	}	
	plots$p_cgc_high <- high_cgc_meth %>% plot_cg_cont_trend(x='`(0.08,0.09]`', y='`(0.09,0.1]`')  + guides(color=FALSE) + coord_cartesian(xlim=c(0,0.3), ylim=c(0,0.3))
	
	plots$p_cgc_chh <- plot_cg_cont_trend(cgc_meth, x='`(0.01,0.02]`', y='CHH', ylab='%CHH') + scale_y_continuous(labels = scales::percent) + guides(color=FALSE) + coord_cartesian(xlim=c(0,1))

	if (length(groups(stats)) > 0){
		legend <- cowplot::get_legend(plots$p_cgc_low)
		plots$p_cgc_low <- plots$p_cgc_low + guides(color=FALSE)
	}	

	if (return_all_figs){
		return(plots)
	}

	p <- cowplot::plot_grid(plotlist=plots, align='h', ncol=NULL, labels="AUTO", label_size=8)

	if (length(groups(stats)) > 0){
		p <- cowplot::plot_grid(p, legend, rel_widths = c(3, .3))	
	}
	
	if (!is.null(ofn)){
		cowplot::save_plot(ofn, p, base_height=base_height, base_aspect_ratio=1, ...)
	}
	
	return(p)
}


# text_summary <- function(db){
# 	n_cells <- nrow(db@cells)
# 	n_cpgs <- nrow(db@cpgs)
# }

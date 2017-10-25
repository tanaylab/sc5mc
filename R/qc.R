#' Plot QC statistics of smat object
#' 
#' 
#' @param cpg_pairs_min_cells minimal number of CpGs to consider in cg_pairs plot
#' @param regions intervals set of capture regions. (default: NULL)
#' @param capture_stats show only statistics from capture regions
#' @param ofn output file
#' @param ... additional parameters for cowplot::save_plot
#' @inheritParams cowplot::save_plot
#' @inheritParams smat.filter_by_cov
#'
#' @return plot
#' 
#' @export
sc5mc.qc_plot <- function(smat, min_cpgs=1, max_cpgs=Inf, min_cells=1, max_cells=Inf, cpg_pairs_min_cells=20, regions=NULL, capture_stats = FALSE, ofn=NULL, width=NULL, base_height=9, ...){
	old <- theme_set(theme_bw(base_size=8))
	on.exit(theme_set(old))
	if (any(min_cpgs != 1, max_cpgs != Inf,  min_cells != 1, max_cells != Inf)){			
		smat <- smat.filter_by_cov(smat, min_cells=min_cells, max_cells=max_cells, max_cpgs=max_cpgs, min_cpgs=min_cpgs)
	}

	if (!is.null(regions) && capture_stats){
		smat <- smat.filter(smat, intervs = regions)
	}

	if (has_stats(smat) && has_name(smat$stats, 'cell_id')){
		smat$stats <- smat$stats %>% mutate(lib = cell_id)
	}

	message('calculating CpG marginals')
	cpg_mars <- sc5mc.plot_cpg_marginals(smat, type='vals') 
	cg_pairs_mars <- sc5mc.plot_cpg_marginals_bars(smat)

	message('calculating cell marginals')
	cell_mars <- sc5mc.plot_cell_marginals(smat, type='vals')

	message('calculating cell pairs')
	cell_pairs_mars <- sc5mc.plot_cell_pairs_coverage(smat)

	figs <- list(cpg_mars, cg_pairs_mars, cell_mars, cell_pairs_mars)	
	
	if (has_stats(smat)){		
		figs[['reads_per_cpg']] <- sc5mc.plot_reads_per_cpg(smat)
		figs[['conversion']] <- sc5mc.plot_conversion(smat)				
	}

	if (has_name(smat, 'tidy_cpgs') || has_name(smat, 'tidy_cpgs_dir')){
		if (!has_name(smat, 'tidy_cpgs')){
			smat <- smat.get_tidy_cpgs(smat)
		}		
		figs[['cells_per_umi']] <- sc5mc.plot_cells_per_umi(smat)
		figs[['reads_per_umi']] <- sc5mc.plot_reads_per_umi(smat, intervals=regions)
		figs[['insert_length']] <- sc5mc.plot_insert_length_reads_per_umi(smat, intervals=regions)

		if (!is.null(regions)){
			figs[['frac_ontar']] <- sc5mc.plot_frac_on_target(smat, intervals=regions)
			figs[['reads_per_region']] <- sc5mc.plot_reads_per_region(smat, intervals=regions)
		}

		# column_num <- 3		
	}
	# } else {
		# column_num <- 2
	# }

	# figs <- map(figs, ~ .x + theme_bw(base_size=8))	
	# p <- cowplot::plot_grid(plotlist=figs, align='none', labels="AUTO", ncol=column_num, label_size=8, vjust=0)
	p <- cowplot::plot_grid(plotlist=figs, align='h', ncol=NULL, labels="AUTO", label_size=8)
	
	if (has_name(smat, 'name')){
		if (capture_stats){
			lab <- glue('{smat$name} (capture)')
		} else {
			lab <- smat$name
		}
		p_title <- cowplot::ggdraw() + cowplot::draw_label(lab, fontface='bold')
		p <- plot_grid(p_title, p, ncol=1, rel_heights=c(0.05, 1, 0.3))
	}	

	if (!is.null(ofn)){
		cowplot::save_plot(ofn, p, base_height=base_height, ...)
	}
	return(p)
}

#' @export
sc5mc.reads_per_cpg <- function(smat){
	if (!has_stats(smat)){
		stop('no stats field!')
	}
	stats <- smat$stats
	if (!has_name(stats, 'cpg_num')){
		cell_mars <- smat.cell_marginals(smat) %>% rename(cpg_num=cov)
		if (has_name(stats, 'lib')){
			stats <- stats %>% rename(cell_id = lib)
		}
		stats <- stats %>% left_join(cell_mars, by='cell_id')
	}	
	return(stats)
}

#' @export
sc5mc.plot_reads_per_cpg <- function(smat){
	stats <- sc5mc.reads_per_cpg(smat)
	
	p <- stats %>% ggplot(aes(x=total_reads, y=cpg_num)) + geom_point(size=0.5, shape=19) + scale_y_continuous(label=comify) + scale_x_continuous(label=scales::comma) + xlab('# of reads') + ylab('# of CpGs') + stat_smooth(method='lm', se=F, linetype='dashed') + ggpmisc::stat_poly_eq(formula = y~ x, eq.with.lhs = "italic(hat(y))~`=`~", aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), size=2, color='darkred', parse = TRUE) + labs(subtitle = qq('median reads = @{comify(round(median(stats$total_reads, na.rm=TRUE)))}\nmedian cpgs = @{comify(round(median(stats$cpg_num, na.rm=TRUE)))}'))
	return(p)	
}

#' @export
sc5mc.plot_conversion <- function(smat){
	if (!has_stats(smat)){
		stop('no stats field!')
	}
	stats <- smat$stats
	p <- stats %>% ggplot(aes(x=CHH)) + geom_density() + scale_x_continuous(labels=scales::percent) + xlab('%C not in CpG context (CHH)')
	return(p)
}


sc5mc.plot_cpg_marginals_bars <- function(smat, cell_nums=c(2,5,10,20,30,50,80)){
	mars <- smat.cpg_marginals(smat)
	purrr::map_dfr(cell_nums, ~ tibble(k=.x, cpgs=sum(mars$cells >= .x))) %>%
		ggplot(aes(x=factor(k), y=log10(cpgs))) + geom_col(width=0.7, fill='darkblue') + scale_y_continuous(labels=comify) + ylab('log10(CpGs with # cells >= x)') + xlab('cells') + labs(subtitle = qq('number of cells = @{comify(ncol(smat$cov))}\nnumber of cpgs = @{comify(nrow(smat$cov))}'))
}


#' plot cpg marginals
#' 
#' @param smat smat object
#'
#' @param type 'percent' or 'vals'
#'
#' @export
sc5mc.plot_cpg_marginals <- function(smat, type='percent'){
	mars <- smat.cpg_marginals(smat) %>% filter(cells > 0) %>% select(cells)
	mars <- mars %>% group_by(cells) %>% summarise(cpgs = n()) %>% arrange(-cells) %>% mutate(c_cpgs = cumsum(cpgs), p = c_cpgs / sum(cpgs))
	x_breaks <- c(1,2,3,4,10,as.numeric(round(quantile(mars$cells, probs=c(0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9, 1))))) %>% unique %>% sort
	if (type == 'percent'){
		p <- mars %>% ggplot(aes(x=cells, y=p)) + geom_point(size=0.3) + scale_x_log10(breaks=x_breaks) + scale_y_continuous(label=scales::percent) + ylab('% of CpGs # cells >= x') + xlab('log10(# of cells)')
	} else {
		p <- mars %>% ggplot(aes(x=cells, y=c_cpgs)) + geom_point(size=0.3) + scale_x_log10(breaks=x_breaks) + scale_y_log10(label=comify) + ylab('log10(CpGs with # cells >= x)') + xlab('log10(# of cells)')
	}
	p <- p + labs(subtitle = qq('number of cells = @{comify(ncol(smat$cov))}\nnumber of cpgs = @{comify(nrow(smat$cov))}')) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
	return(p)
}

#' plot cell marginals
#' 
#' @param smat smat object
#'
#' @param type 'percent' or 'vals'
#'
#' @export
sc5mc.plot_cell_marginals <- function(smat, type='percent'){
	mars <- smat.cell_marginals(smat) %>% filter(cov > 0) %>% select(cov)
	mars <- mars %>% group_by(cov) %>% summarise(cells = n())  %>% arrange(-cov) %>% mutate(c_cells = cumsum(cells), p = c_cells / sum(cells))

	x_breaks <- as.numeric(round(quantile(mars$cov, probs=c(0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9, 1))))

	if (type == 'percent'){
		p <-mars %>% ggplot(aes(x=cov, y=p)) + geom_point(size=0.3) + scale_x_log10(label=comify, breaks=x_breaks) + scale_y_continuous(label=scales::percent) + ylab('% cells # CpGs >= x') + xlab('log10(# of CpGs)')
	} else {
		p <- mars %>% ggplot(aes(x=cov, y=c_cells)) + geom_point(size=0.3) + scale_x_log10(label=comify) + scale_y_continuous(label=comify) + ylab('cells with # CpGs >= x') + xlab('log10(# of CpGs)')
	}

	p <- p + labs(subtitle = qq('number of cells = @{comify(ncol(smat$cov))}\nnumber of cpgs = @{comify(nrow(smat$cov))}'))
	return(p)
}

#' plot cell pairs coverage
#' 
#' @param smat smat object
#'
#' @export
sc5mc.plot_cell_pairs_coverage <- function(smat){
	mars <- smat.cell_pairs_marginals(smat)
	cmars <- mars %>% count(ntot) %>% arrange(-ntot) %>% mutate(cn = cumsum(n), p = cn / sum(n))

	x_breaks <- as.numeric(round(quantile(cmars$ntot, probs=c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1))))
	p <- cmars %>% ggplot(aes(x=ntot, y=p)) + geom_point(size=0.3) + xlab('log10(# of jointly covered CpGs)') + ylab('% of cell pairs') + scale_x_log10(label=comify, breaks=x_breaks) + scale_y_continuous(labels=scales::percent) + labs(subtitle = qq('number of cells = @{comify(length(unique(mars$cell1)))}\nnumber of pairs = @{comify(nrow(mars))}')) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
	return(p)
}

#' plot cg pairs coverage
#' 
#' @param smat smat object
#'
#' @export
sc5mc.plot_cg_pairs_coverage <- function(smat, min_cells=30){
	mars <- smat.cpg_pairs_marginals(smat,  min_cells=min_cells)
	cmars <- mars %>% arrange(-ncells) %>% mutate(cn = cumsum(npairs), p = cn / sum(npairs))

	ncpgs <- sum(rowSums(smat$cov) >= min_cells)

	x_breaks <- as.numeric(round(quantile(cmars$ncells, probs=c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1))))
	
	p <- cmars %>% ggplot(aes(x=ncells, y=p)) + geom_point(size=0.3) + xlab('log10(# of jointly covered cells)') + ylab('% of CpG pairs') + scale_x_log10(label=comify, breaks=x_breaks) + scale_y_continuous(labels=scales::percent) + labs(subtitle = qq('number of cpgs = @{comify(ncpgs)} (covered by at least @{min_cells} cells)\nnumber of pairs = @{comify(sum(mars$npairs))}'))
	
	return(p)
}

#' plot indexes statistics
#' 
#' @param smat smat object
#' @param ofn figure filename (ggsave parameter)
#' @param width figure width (ggsave parameter)
#' @param height figure height (ggsave parameter)
#' @param cells_per_umi plot cells per umi
#' 
#' @return plots of total reads, %mapping and %unique 
#'
#' @export
sc5mc.index_plots <- function(smat, ofn, width=7*2.5, height=5*2.3, cells_per_umi=FALSE){
	stats <- smat$stats
    stats <- stats %>% mutate(row = factor(row, levels= rev(sort(unique(as.character(row))))), column = factor(column, levels=sort(unique(as.numeric(column)))))

    p_frac_mapped <- stats %>% ggplot(aes(x=column, y=row, color=empty, fill=mapped_frac)) + geom_tile(size=0.5)+ scale_color_manual(values=rev(c('red', 'gray'))) + scale_fill_gradientn(colors=c('white', 'darkgreen')) + guides(color=F) 

    p_tot_reads <- stats %>% ggplot(aes(x=column, y=row, color=empty, fill=log10(total_reads))) + geom_tile(size=0.5) + scale_color_manual(values=rev(c('black', 'gray'))) + scale_fill_gradientn(colors=c('white', 'white', 'red')) + guides(color=F) 

    p_uniq_frac <- stats %>% ggplot(aes(x=column, y=row, color=empty, fill=uniq_frac)) + geom_tile(size=0.5) + scale_color_manual(values=rev(c('black', 'gray'))) + scale_fill_gradientn(colors=c('white', 'red')) + guides(color=F) 

    if (cells_per_umi){
    	p_cells_per_umi <- sc5mc.plot_cells_per_umi(smat)	
    	p <- cowplot::plot_grid(p_tot_reads, p_uniq_frac, p_frac_mapped, p_cells_per_umi, align='hv', ncol=2) 
    } else {
    	p <- cowplot::plot_grid(p_tot_reads, p_uniq_frac, p_frac_mapped, align='hv', ncol=2) 
    }    

    p + ggsave(ofn, width=width, height=height)
}

#' Calculates cells per UMI
#' 
#' @param smat smat object
#' 
#' @return tibble with number of cells ('n_cells') and number of molecules ('n')
#'
#' @export
sc5mc.cells_per_umi <- function(smat){
	if (!has_name(smat, 'tidy_cpgs_all')){
		message("Getting tidy cpgs")
		smat <- smat.get_tidy_cpgs(smat, unique=FALSE)
	}
	message("Calculating reads per UMI")
	tcpgs <- smat$tidy_cpgs  %>% distinct(read_id, .keep_all=TRUE) %>% select(cell_id, read_id, chrom, start, end, strand, umi1, umi2, insert_len)
	cpu <- tcpgs %>% filter(end != '-') %>% group_by(chrom, start, end) %>% summarise(n = n(), n_cells=n_distinct(cell_id)) %>% ungroup()
	cpu <- cpu %>% group_by(n_cells) %>% summarise(n = n()) %>% mutate(p = n / sum(n))
	return(cpu)	 
}

#' Plots cells per UMI
#' 
#' @param smat smat object
#' 
#' @return plot of number of cells versus percent of moleculs (only for molecules that appeared in more than 1 cell)
#'
#' @export
sc5mc.plot_cells_per_umi <- function(smat){
	cpu <- sc5mc.cells_per_umi(smat)	
	p <- cpu %>% filter(n_cells > 1) %>% ggplot(aes(x=factor(n_cells), y=p)) + geom_col(fill='darkblue') + scale_y_continuous(labels=scales::percent) + xlab('# of cells') + ylab('% of molecules') + labs(subtitle='') + coord_cartesian(ylim=c(0,0.01))
	return(p)
}

#' @export
sc5mc.ontar_reads <- function(smat, intervals=NULL, return_orig=TRUE){
	if (!has_name(smat, 'tidy_cpgs')){
		message("Getting tidy cpgs")
		smat <- smat.get_tidy_cpgs(smat, unique=TRUE)
	}
	
	if (!is.null(intervals)){
		reads <- smat$tidy_cpgs %>% distinct(read_id, num, .keep_all=T) %>% filter(end != '-')
		tcpgs <- .gpatterns.filter_ontar_reads(reads, intervals, return_orig=return_orig)		
	} else {
		tcpgs <- smat$tidy_cpgs
	}
	return(tcpgs)
}

#' @export
sc5mc.reads_per_umi <- function(smat, intervals=NULL){
	tcpgs <- sc5mc.ontar_reads(smat, intervals)	
	tcpgs <- tcpgs  %>% distinct(read_id, num) 

	rpu <- tcpgs %>% group_by(reads=num) %>% summarise(n=n())
	return(rpu)	
}

#' @export
sc5mc.plot_reads_per_umi <- function(smat, intervals=NULL){
	rpu <- sc5mc.reads_per_umi(smat, intervals)
	rpu <- rpu %>% 
		mutate(reads = cut(reads, breaks=c(0,1,4,10,1e6), include.lowest=TRUE, labels=c('1', '2-4', '5-10', '>10'))) %>% 
		group_by(reads) %>% 
		summarise(n = sum(n, na.rm=TRUE)) %>% 
		mutate(p = n / sum(n))	
	p <- rpu %>% ggplot(aes(x=reads, y=p)) + geom_col(width=0.7, fill='midnightblue') + scale_y_continuous(label=scales::percent) +  xlab('Reads / UMI') + ylab("% of reads") +theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5)) + labs(subtitle=glue("number of reads: {comify(sum(rpu$n))}"))
	return(p)	
}

#' @export
sc5mc.plot_insert_length_reads_per_umi <- function(smat, intervals=NULL){
	tcpgs <- sc5mc.ontar_reads(smat, intervals)	
 	rpu <- tcpgs %>% mutate(reads = cut(num, breaks=c(0,1,4,10,1e6), include.lowest=TRUE, labels=c('1', '2-4', '5-10', '>10'))) %>% filter(!is.na(reads))
 	p <- rpu %>% ggplot(aes(x=abs(insert_len), color=reads)) + geom_density(size=0.9) + xlab('Insert length') +ylab('Density') + ggsci::scale_color_aaas(name = 'Reads / UMI') + theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5), legend.position = c(0.8, 0.8), legend.box.background = element_rect()) + coord_cartesian(xlim=c(0,450)) + guides(color=guide_legend(keywidth=0.5, keyheight=0.5))
	
	return(p)	
}


#' @export
sc5mc.frac_on_target <- function(smat, intervals){
	if (!has_name(smat, 'tidy_cpgs')){
		message("Getting tidy cpgs")
		smat <- smat.get_tidy_cpgs(smat, unique=TRUE)
	}

	reads <- smat$tidy_cpgs %>% distinct(read_id, .keep_all=T) %>% filter(end != '-')
	ontar <- .gpatterns.filter_ontar_reads(reads, intervals)	
	stats <- reads %>% group_by(cell_id) %>% summarise(all_reads = sum(num)) %>% left_join(ontar %>%  group_by(cell_id) %>% summarise(ontar_reads = sum(num)), by='cell_id') %>% replace_na(ontar_reads = 0) %>% mutate(frac_ontar = ontar_reads / all_reads)
	
	return(stats)	
}

#' @export
sc5mc.plot_frac_on_target <- function(smat, intervals, min_reads=1e3){
	stats <- sc5mc.frac_on_target(smat, intervals)	
	p <- stats %>% filter(all_reads >= min_reads) %>% ggplot(aes(x=frac_ontar)) + geom_histogram(binwidth=0.01, fill='darkblue') + ylab("# of cells") + xlab("% reads 'on target'") + scale_x_continuous(labels=scales::percent) + scale_y_continuous(labels=comify) + labs(subtitle=glue("Overall % reads on target: {scales::percent(sum(stats$ontar_reads, na.rm=TRUE) / sum(stats$all_reads, na.rm=TRUE))}"))
	return(p)	
}

#' @export
sc5mc.reads_per_region <- function(smat, intervals){
	reads <- sc5mc.ontar_reads(smat, intervals, return_orig=FALSE)	
   	covs <- reads %>% 
        gintervals.neighbors1(intervals) %>% 
        group_by(chrom1, start1, end1) %>% 
        summarise(n = n()) %>% 
        select(chrom=chrom1, start=start1, end=end1, cov=n)

    covs <- covs %>%
        full_join(intervals, by=c('chrom', 'start', 'end')) %>% 
        replace_na(replace=list(cov = 0)) %>% 
		ungroup()

    return(covs) 	
}

#' @export
sc5mc.plot_reads_per_region <- function(smat, intervals){
	reg_stats <- sc5mc.reads_per_region(smat, intervals)	
	p <- reg_stats %>% 
		ungroup() %>% 
		count(cov) %>% 
		mutate(cum = rev(cumsum(rev(n)))) %>% 
		ggplot(aes(x=cov, y=cum)) + 
			geom_point(size=0.5, shape=19) + 
			xlab("Reads") + 
			ylab("Regions with at least x reads") + 
			scale_y_continuous(labels=comify) +
			labs(subtitle=glue("number of reads: {comify(sum(reg_stats$cov))}"))
	return(p) 
}


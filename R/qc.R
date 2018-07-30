#' Plot QC statistics of smat object
#' 
#' @param smat smat object
#' @param cpg_pairs_min_cells minimal number of CpGs to consider in cg_pairs plot
#' @param regions intervals set of capture regions. (default: NULL)
#' @param capture_stats show only statistics from capture regions
#' @param ofn output file
#' @param subtitle additional text for title
#' @param raw_reads_dir directory with the raw reads
#' @param db cgdb object (in order to plot CpG conent / average methylation plots)
#' @param min_cgc_cov minimal coverage for average methylation per CpG content estimation
#' @param ... additional parameters for cowplot::save_plot
#' @inheritParams cowplot::save_plot
#' @inheritParams smat.filter_by_cov
#'
#' @return plot
#' 
#' @export
sc5mc.qc_plot <- function(smat, min_cpgs=1, max_cpgs=Inf, min_cells=1, max_cells=Inf, cpg_pairs_min_cells=20, regions=NULL, capture_stats = FALSE, ofn=NULL, width=NULL, base_height=9, subtitle = NULL, raw_reads_dir=NULL, db=NULL, min_cgc_cov=300, ...){
	opt <- options(gmax.data.size=1e9)
	on.exit(options(opt))

	old <- theme_set(theme_bw(base_size=8))
	on.exit(theme_set(old))
	if (any(min_cpgs != 1, max_cpgs != Inf,  min_cells != 1, max_cells != Inf)){			
		smat <- smat.filter_by_cov(smat, min_cells=min_cells, max_cells=max_cells, max_cpgs=max_cpgs, min_cpgs=min_cpgs)
	}

	if (!is.null(regions)){
		smat_ontar <- smat.filter(smat, intervs = regions)
		if (capture_stats){
			smat <- smat_ontar
		}		
	}

	if (has_stats(smat) && has_name(smat$stats, 'cell_id')){
		smat$stats <- smat$stats %>% mutate(lib = cell_id)
	}

	batch_colors <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')

	index_plots <- sc5mc.index_plots(smat, return_plotlist = TRUE)
	index_plots <- map(index_plots, ~ .x + theme(legend.position="bottom"))	

	figs <- index_plots

	# stats_fn <- glue('{smat$path}_stats.csv')
	if (has_stats(smat) && !has_name(smat$stats, 'illu_read_num') && !is.null(raw_reads_dir)){
		smat$stats <- smat.raw_reads_stats(smat, raw_reads_dir)
		all_reads <- smat$stats %>% distinct(illumina_index, illu_read_num) %>% pull(illu_read_num) %>% sum(na.rm=TRUE)
		figs[['p_demulti']] <- smat$stats %>%
			group_by(batch_id) %>% 
			summarise(total_reads = sum(total_reads, na.rm=TRUE)) %>% 
			bind_rows(smat$stats  %>% 
				summarise(batch_id = 'no_index', total_reads = all_reads - sum(total_reads, na.rm=TRUE))) %>% 
			ungroup() %>% 
			mutate(p = total_reads / sum(total_reads, na.rm=TRUE)) %>% 
			ggplot(aes(x=batch_id, y= p, fill=batch_id)) + 
				geom_col() + 
				scale_fill_manual(values=c(batch_colors, 'black'), drop=FALSE) + 
				guides(fill=FALSE) + 
				xlab('') + 
				ylab('% of reads') + 
				scale_y_continuous(label=scales::percent) + 
				labs(subtitle = glue('total number of reads: {scales::comma(all_reads)}'))
		
		# fwrite(stats, stats_fn)
	}

	
	if (has_stats(smat)){	
		loginfo('calculating reads per CpG')
		
		figs[['reads_per_cpg']] <- sc5mc.plot_reads_per_cpg(smat, color_column='batch_id', colors=batch_colors, log_scale=TRUE) + guides(color=FALSE)
	}

	loginfo('calculating CpG marginals')
	figs[['cpg_mars']] <- sc5mc.plot_cpg_marginals(smat, type='vals') 
	figs[['cg_pairs_bars']] <- sc5mc.plot_cpg_marginals_bars(smat)

	loginfo('calculating cell marginals')
	figs[['cell_mars']] <- sc5mc.plot_cell_marginals(smat, type='vals')
	figs[['cell_mars_bars']] <- sc5mc.plot_cell_marginals_bars(smat)

	if (!is.null(regions)){
		figs[['cell_mars_bars_ontar']] <- sc5mc.plot_cell_marginals_bars(smat_ontar) + xlab('On target CpGs')
	}

	# loginfo('calculating cell pairs')
	# cell_pairs_mars <- sc5mc.plot_cell_pairs_coverage(smat)	
	
	if (has_stats(smat)){			
		loginfo('calculating conversion')
		figs[['conversion']] <- sc5mc.plot_conversion(smat)				
	}

	if (has_name(smat, 'tidy_cpgs') || has_name(smat, 'tidy_cpgs_dir')){

		smat.cache_stats(smat, regions=regions)

		glob_stats_regions <- NULL
		if (capture_stats){
			glob_stats_regions <- regions
		} 

		if (!is.null(regions)){
			loginfo('calculating reads per umi')		
			figs[['reads_per_umi']] <- sc5mc.plot_capture_reads_per_umi(smat, intervals=regions)			
		} else {
			loginfo('calculating reads per umi')		
			figs[['reads_per_umi']] <- sc5mc.plot_reads_per_umi(smat, intervals=glob_stats_regions)	
		}		
	
		loginfo('calculating reads per umi per insert length')
		figs[['insert_length']] <- sc5mc.plot_insert_length_reads_per_umi(smat, intervals=glob_stats_regions)
	
		loginfo('calculating cells per umi')
		figs[['cells_per_umi']] <- sc5mc.plot_cells_per_umi(smat)
		
		if (!is.null(regions)){
			loginfo('calculating on target')
			figs[['frac_ontar']] <- sc5mc.plot_frac_on_target(smat, intervals=regions)
			loginfo('calculating reads per region')
			figs[['reads_per_region']] <- sc5mc.plot_reads_per_region(smat, intervals=regions)
		}

	}

	if (!is.null(db)){
		db_figs <- qc_plot(db, min_cgc_cov=min_cgc_cov, return_all_figs=TRUE)
		figs[['cgc_low']] <- db_figs[['p_cgc_low']]
		figs[['cgc_high']] <- db_figs[['p_cgc_high']]
		figs[['cgc_low_high']] <- db_figs[['p_cgc_low_high']]		
	}
	
	loginfo('plotting...')
	png(tempfile(fileext='.png'))
	p <- cowplot::plot_grid(plotlist=figs, align='h', ncol=NULL, labels="AUTO", label_size=8)
	dev.off()	
	
	if (has_name(smat, 'name')){
		if (capture_stats){
			lab <- glue('{smat$name} (capture)')
		} else {
			lab <- smat$name
		}
		
		if (!is.null(subtitle)){
			lab <- paste(lab, subtitle, sep='\n')
		}
	
		p_title <- cowplot::ggdraw() + cowplot::draw_label(lab, fontface='bold')	
		p <- cowplot::plot_grid(p_title, p, ncol=1, rel_heights=c(0.05, 1, 0.3))
	}	

	if (!is.null(ofn)){
		loginfo('saving plot')
		cowplot::save_plot(ofn, p, base_height=base_height, ...)
	}

	return(p)
}

#' @export
smat.cache_stats <- function(smat, regions=NULL, overwrite=FALSE){
	cpu_fn <- glue('{smat$tidy_cpgs_dir}/stats/cpu.csv')
	rpu_fn <- glue('{smat$tidy_cpgs_dir}/stats/rpu.csv')
	rpu_ontar_fn <- glue('{smat$tidy_cpgs_dir}/stats/rpu_ontar.csv')
	rpu_offtar_fn <- glue('{smat$tidy_cpgs_dir}/stats/rpu_offtar.csv')
	rpu_insert_length_fn <- glue('{smat$tidy_cpgs_dir}/stats/rpu_insert_length.csv')
	ontar_stats_fn <- glue('{smat$tidy_cpgs_dir}/stats/ontar_stats.csv')
	reg_stats_fn <- glue('{smat$tidy_cpgs_dir}/stats/regions_stats.csv')


	if (!overwrite && all(file.exists(c(cpu_fn, rpu_fn, rpu_insert_length_fn)))){
		if (is.null(regions) || all(file.exists(c(ontar_stats_fn, reg_stats_fn, rpu_ontar_fn, rpu_offtar_fn)))){
			return(invisible(NULL))
		}		
	}	

	loginfo('caching statistics')

	if (!has_name(smat, 'tidy_cpgs')){
		loginfo("Getting tidy cpgs unique")
		smat <- smat.get_tidy_cpgs(smat, unique=TRUE)
	}		
	if (!has_name(smat, 'tidy_cpgs_all')){
		loginfo("Getting tidy cpgs non unique")
		smat <- smat.get_tidy_cpgs(smat, unique=FALSE)
	}

	dir.create(glue('{smat$tidy_cpgs_dir}/stats'), showWarnings=FALSE)

	if (!file.exists(cpu_fn) || overwrite){
		cpu <- sc5mc.cells_per_umi(smat)
		fwrite(cpu, cpu_fn)	
	}
		
	if (!file.exists(rpu_fn) || overwrite){
		rpu <- sc5mc.reads_per_umi(smat)	
		fwrite(rpu, rpu_fn)
	}
	

	if (!file.exists(rpu_insert_length_fn) || overwrite){
		rpu_insert_length <- sc5mc.reads_per_umi_insert_length(smat)
		fwrite(rpu_insert_length, rpu_insert_length_fn)
	}

	if (!is.null(regions)){
		
		if (!file.exists(ontar_stats_fn) || overwrite){
			ontar_stats <- sc5mc.frac_on_target(smat, regions)			
			fwrite(ontar_stats, ontar_stats_fn)
		}

		if (!file.exists(reg_stats_fn) || overwrite){
			reg_stats <- sc5mc.reads_per_region(smat, regions)	
			fwrite(reg_stats, reg_stats_fn)	
		}

		if (!file.exists(rpu_ontar_fn) || overwrite){
			rpu_ontar <- sc5mc.reads_per_umi(smat, intervals=regions)	
			fwrite(rpu_ontar, rpu_ontar_fn)
		}

		if (!file.exists(rpu_offtar_fn) || overwrite){
			rpu_offtar <- sc5mc.reads_per_umi(smat, intervals=gintervals.diff(gintervals.all(), regions))
			fwrite(rpu_offtar, rpu_offtar_fn)
		}		
	}
}

smat.raw_reads_stats <- function(smat, raw_reads_dir){
	stats <- smat$stats
	
	raw_fns <- map_df(unique(stats$illumina_index), ~ data.frame(illumina_index = .x) %>% mutate(fn = list(list.files(file.path(raw_reads_dir, .x, 'raw'), full.names=TRUE, pattern='R1.*\\.fastq\\.gz'))) ) %>% unnest(fn)
		
	read_stats <- plyr::adply(raw_fns, 1, function(.x) {		
		.x %>% mutate(read_num = system(glue('zcat {.x$fn} | sed -n 2~4p | wc -l'), intern=TRUE))
	}, .parallel = TRUE)

	read_stats <- read_stats %>% group_by(illumina_index) %>% summarise(illu_read_num = sum(as.numeric(read_num), na.rm=TRUE) )
	
	stats <- stats %>% left_join(read_stats, by = "illumina_index")

	return(stats)
}

# count_dir_reads <- function(directory, illumina_index){
# 	fastq_files <- list.files(file.path(directory, illumina_index, 'raw'), full.names=TRUE, pattern='R1.*\\.fastq\\.gz')
# 	browser()

# 	read_num <- system(glue('zcat {paste(fastq_files, collapse=" ")} | sed -n 2~4p | wc -l'), intern=TRUE)

# }

#' @export
sc5mc.reads_per_cpg <- function(smat){
	if (!has_stats(smat)){
		stop('no stats field!')
	}
	stats <- smat$stats
	if (!has_name(stats, 'cpg_num')){
		cell_mars <- smat.cell_marginals(smat) %>% rename(cpg_num=cov)
		if (has_name(stats, 'lib')  && !has_name(stats, 'cell_id')){
			stats <- stats %>% rename(cell_id = lib)
		}
		stats <- stats %>% left_join(cell_mars, by='cell_id')
	}	
	return(stats)
}

#' @export
sc5mc.plot_reads_per_cpg <- function(smat, color_column = NULL, colors = NULL, log_scale=FALSE){
	if (class(smat) == 'smat'){
		stats <- sc5mc.reads_per_cpg(smat)	
	}

	if (class(smat) == 'cgdb'){
		stats <- smat@cells
		if (length(groups(stats)) > 0){
			color_column <- color_column %||% 'type'
			stats <- stats %>% unite_(color_column, groups(stats))
		}
		if (has_name(stats, 'cg_num')){
			stats$cpg_num <- stats$cg_num
		}
	}
	max_read_num <- max(stats$total_reads, na.rm=TRUE)

	log_lines_d <- bind_rows(tibble(x = c(1,max_read_num)) %>% mutate(y=0.1*x, model='0.1'), 
                             tibble(x = c(1,max_read_num)) %>% mutate(y=0.2*x, model='0.2'), 
                             tibble(x = c(1,max_read_num)) %>% mutate(y=0.3*x, model='0.3'),
                             tibble(x = c(1,max_read_num)) %>% mutate(y=0.4*x, model='0.4'))

	p <- stats %>% 
		ggplot(aes_string(x='total_reads', y='cpg_num', color=color_column, shape='empty', size='empty')) + 
			# geom_point(size=0.5, shape=19) + 
			geom_point() + 
			scale_shape_manual(values=c(19,4)) +
			 scale_size_manual(values=c(0.5, 2)) +
			# stat_smooth(method='lm', se=F, linetype='dashed') + 
			guides(shape=FALSE, size=FALSE) + 
			geom_line(d = log_lines_d, inherit.aes = FALSE, aes(x=x, y=y, group=model), linetype='dashed', color='darkgray') +
			labs(subtitle = qq('median reads = @{comify(round(median(stats$total_reads[stats$cpg_num >= 1000], na.rm=TRUE)))}\nmedian cpgs = @{comify(round(median(stats$cpg_num[stats$cpg_num >= 1000], na.rm=TRUE)))}'))

	if (log_scale){
		p <- p + scale_y_log10(label=scales::scientific) + scale_x_log10(label=scales::scientific) + annotation_logticks()
	} else {
		p <- p + scale_y_continuous(label=scales::scientific) + scale_x_continuous(label=scales::scientific) + xlab('# of reads') + ylab('# of CpGs')
	}	
	
	if (!is.null(color_column)){
		if (is.null(colors)){
			p <- p + ggsci::scale_color_lancet()	
		} else {
			p <- p + scale_color_manual(values = colors)
		}
		
	} else {
		# p <- p + ggpmisc::stat_poly_eq(formula = y~ x, eq.with.lhs = "italic(hat(y))~`=`~", aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), size=2, color='darkred', parse = TRUE)
	}
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
	if (class(smat) == 'smat'){
		mars <- smat.cpg_marginals(smat)
		n_cells <- ncol(smat$cov)
		n_cpgs <- nrow(smat$cov)
	}

	if (class(smat) == 'cgdb'){
		mars <- smat %>% ungroup_cells() %>% ungroup_cpgs() %>% summarise_cells() %>% filter(cov > 0) %>% select(cells=cov)
		n_cells <- nrow(smat@cells)
		n_cpgs <- nrow(smat@cpgs)
	}		

	p <- purrr::map_dfr(cell_nums, ~ tibble(k=.x, cpgs=sum(mars$cells >= .x))) %>%
		ggplot(aes(x=factor(k), y=cpgs)) + geom_col(width=0.7, fill='darkblue') + ylab('CpGs with # cells >= x') + scale_y_log10(labels=comify) + xlab('cells') +labs(subtitle = qq('number of cells = @{comify(n_cells)}\nnumber of cpgs = @{comify(n_cpgs)}'))
	return(p)
}


#' plot cpg marginals
#' 
#' @param smat smat object / cgdb object
#'
#' @param type 'percent' or 'vals'
#'
#' @export
sc5mc.plot_cpg_marginals <- function(smat, type='percent'){
	if (class(smat) == 'smat'){
		mars <- smat.cpg_marginals(smat) %>% filter(cells > 0) %>% select(cells)
		n_cells <- ncol(smat$cov)
		n_cpgs <- nrow(smat$cov)
	}

	if (class(smat) == 'cgdb'){
		mars <- smat %>% summarise_cells() %>% filter(cov > 0) %>% select(cells=cov)
		n_cells <- nrow(smat@cells)
		n_cpgs <- nrow(smat@cpgs)
	}	
	
	mars <- mars %>% group_by(cells) %>% summarise(cpgs = n()) %>% arrange(-cells) %>% mutate(c_cpgs = cumsum(cpgs), p = c_cpgs / sum(cpgs))
	x_breaks <- c(1,2,3,4,10,as.numeric(round(quantile(mars$cells, probs=c(0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9, 1))))) %>% unique %>% sort
	if (type == 'percent'){
		p <- mars %>% ggplot(aes(x=cells, y=p)) + geom_point(size=0.3) + scale_x_log10(breaks=x_breaks) + scale_y_continuous(label=scales::percent) + ylab('% of CpGs # cells >= x') + xlab('log10(# of cells)')
	} else {
		p <- mars %>% ggplot(aes(x=cells, y=c_cpgs)) + geom_point(size=0.3) + scale_x_log10(breaks=x_breaks) + scale_y_log10(label=comify) + ylab('log10(CpGs with # cells >= x)') + xlab('log10(# of cells)')
	}
	p <- p + labs(subtitle = qq('number of cells = @{comify(n_cells)}\nnumber of cpgs = @{comify(n_cpgs)}')) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
	return(p)
}

sc5mc.plot_cell_marginals_bars <- function(smat, cpg_nums = c(500, 1000, 5000, 10000, 5e4, 8e4, 1e5, 1.5e5, 2e5)){
	if (class(smat) == 'smat'){
		cell_mars <- smat.cell_marginals(smat) %>% filter(cov > 0) %>% select(cov)
		n_cells <- ncol(smat$cov)
		n_cpgs <- nrow(smat$cov)
	} 

	if (class(smat) == 'cgdb'){
		cell_mars <- smat %>% ungroup_cells() %>% ungroup_cpgs() %>% summarise_cpgs() %>% filter(cov > 0) %>% select(cov)
		n_cells <- nrow(smat@cells)
		n_cpgs <- nrow(smat@cpgs)		
	}

	p <- purrr::map_dfr(cpg_nums, ~ tibble(k=.x, cpgs=sum(cell_mars$cov >= .x))) %>% ggplot(aes(x=factor(k), y=cpgs)) + geom_col(width=0.7, fill='darkblue') + ylab('cells with # CpGs >= x') + xlab('CpGs') + scale_x_discrete(labels=function(x) scales::comma(as.numeric(x))) + labs(subtitle = qq('number of cells = @{comify(n_cells)}\nnumber of cpgs = @{comify(n_cpgs)}')) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
	return(p)
}



#' plot cell marginals
#' 
#' @param smat smat object / cgdb object
#'
#' @param type 'percent' or 'vals'
#'
#' @export
sc5mc.plot_cell_marginals <- function(smat, type='percent'){
	if (class(smat) == 'smat'){
		mars <- smat.cell_marginals(smat) %>% filter(cov > 0) %>% select(cov)
		n_cells <- ncol(smat$cov)
		n_cpgs <- nrow(smat$cov)
	} 

	if (class(smat) == 'cgdb'){
		mars <- smat %>% summarise_cpgs() %>% filter(cov > 0) %>% select(cov)
		n_cells <- nrow(smat@cells)
		n_cpgs <- nrow(smat@cpgs)
	}
	
	mars <- mars %>% group_by(cov) %>% summarise(cells = n())  %>% arrange(-cov) %>% mutate(c_cells = cumsum(cells), p = c_cells / sum(cells))

	x_breaks <- as.numeric(round(quantile(mars$cov, probs=c(0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9, 1))))

	if (type == 'percent'){
		p <-mars %>% ggplot(aes(x=cov, y=p)) + geom_point(size=0.3) + scale_x_log10(label=comify, breaks=x_breaks) + scale_y_continuous(label=scales::percent) + ylab('% cells with # CpGs >= x') + xlab('log10(# of CpGs)')
	} else {
		p <- mars %>% ggplot(aes(x=cov, y=c_cells)) + geom_point(size=0.3) + scale_x_log10(label=comify) + scale_y_continuous(label=comify) + ylab('cells with # CpGs >= x') + xlab('log10(# of CpGs)')
	}

	p <- p + labs(subtitle = qq('number of cells = @{comify(n_cells)}\nnumber of cpgs = @{comify(n_cpgs)}'))
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
#' @param return_plotlist return a list of plots
#' 
#' @return plots of total reads, %mapping and %unique 
#'
#' @export
sc5mc.index_plots <- function(smat, ofn=NULL, width=7*2.5, height=5*2.3, cells_per_umi=FALSE, return_plotlist=FALSE, min_reads=1e3){
	stats <- smat$stats
    stats <- stats %>% mutate(row = factor(row, levels= rev(sort(unique(as.character(row))))), column = factor(column, levels=sort(unique(as.numeric(column)))))
    stats <- stats %>% mutate(batch_id = ifelse(empty, 'empty', batch_id), batch_id = factor(batch_id, levels = c(paste0('batch', 1:4), 'empty')))

    stats <- stats %>% mutate(read_per_cpg = total_reads / cg_num)
    
    # p_frac_mapped <- stats %>% ggplot(aes(x=column, y=row, color=empty, fill=mapped_frac)) + geom_tile(size=0.5)+ scale_color_manual(values=rev(c('red', 'gray'))) + scale_fill_gradientn(colors=c('white', 'darkgreen')) + guides(color=F) 

    p_cpg_num <- stats %>% mutate(cg_num = ifelse(total_reads >= 1e3, cg_num, NA), cg_num = ifelse(cg_num >= 1e6, 1e6, cg_num)) %>% ggplot(aes(x=column, y=row, color=empty, fill=log10(cg_num))) + geom_tile(size=0.5)+ scale_color_manual(values=rev(c('black', 'gray'))) + scale_fill_gradientn(colors=c('white', 'white', 'red', 'yellow'), limits=c(1, 6)) + guides(color=F) 

    p_tot_reads <- stats %>% mutate(total_reads = ifelse(total_reads >= 1e3, total_reads, NA), total_reads = ifelse(total_reads >= 5e6, 5e6, total_reads)) %>% ggplot(aes(x=column, y=row, color=empty, fill=log10(total_reads))) + geom_tile(size=0.5) + scale_color_manual(values=rev(c('black', 'gray'))) + scale_fill_gradientn(colors=c('white', 'white', 'red', 'yellow'), limits=c(3, log10(5e6))) + guides(color=F) 

    p_read_per_cpg <- stats %>% mutate(read_per_cpg = ifelse(total_reads >= 1e3, read_per_cpg, NA)) %>% ggplot(aes(x=column, y=row, color=empty, fill=log2(read_per_cpg))) + geom_tile(size=0.5) + scale_color_manual(values=rev(c('black', 'gray'))) + scale_fill_gradientn(colors=c('white', 'white', 'red')) + guides(color=F) 

    # p_uniq_frac <- stats %>% ggplot(aes(x=column, y=row, color=empty, fill=uniq_frac)) + geom_tile(size=0.5) + scale_color_manual(values=rev(c('black', 'gray'))) + scale_fill_gradientn(colors=c('white', 'red')) + guides(color=F) 

    batch_colors <- c('#e41a1c','#377eb8','#4daf4a','#984ea3', 'black')
    p_batches <- stats %>% ggplot(aes(x=column, y=row, fill=batch_id)) + geom_tile(size=0.5) + scale_fill_manual(values=batch_colors, drop=FALSE)

    
    if (return_plotlist){
    	return(list(p_batches = p_batches, p_tot_reads=p_tot_reads, p_cpg_num=p_cpg_num, p_read_per_cpg=p_read_per_cpg))
    }

    if (cells_per_umi){
    	p_cells_per_umi <- sc5mc.plot_cells_per_umi(smat)	
    	p <- cowplot::plot_grid(p_batches, p_tot_reads, p_cpg_num, p_read_per_cpg, p_cells_per_umi, align='hv', ncol=2) 
    } else {
    	p <- cowplot::plot_grid(p_batches, p_tot_reads, p_cpg_num, p_read_per_cpg, align='hv', ncol=2) 
    }    

    if (!is.null(ofn)){
    	p <- p + ggsave(ofn, width=width, height=height)	
    }

    return(p)

    
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
		loginfo("Getting tidy cpgs")
		smat <- smat.get_tidy_cpgs(smat, unique=FALSE)
	}
	loginfo("Calculating reads per UMI")
	tcpgs <- smat$tidy_cpgs_all %>% distinct(read_id, .keep_all=TRUE) %>% select(cell_id, read_id, chrom, start, end, strand, umi1, umi2, insert_len)
	cpu <- tcpgs %>% filter(end != '-') %>% group_by(chrom, start, end) %>% summarise(n = n(), n_cells=n_distinct(cell_id)) %>% ungroup()
	cpu <- cpu %>% group_by(n_cells) %>% summarise(n = n()) %>% mutate(p = n / sum(n))
	return(cpu)	 
}

#' Plots cells per UMI
#' 
#' @param smat smat object
#' 
#' @return plot of number of cells versus percent of molecules (only for molecules that appeared in more than 1 cell)
#'
#' @export
sc5mc.plot_cells_per_umi <- function(smat){
	cpu_fn <- glue('{smat$tidy_cpgs_dir}/stats/cpu.csv')
	if (file.exists(cpu_fn)){
		cpu <- fread(cpu_fn) %>% as_tibble()
	} else {
		cpu <- sc5mc.cells_per_umi(smat)					
	}	
	p <- cpu %>% filter(n_cells > 1) %>% ggplot(aes(x=factor(n_cells), y=p)) + geom_col(fill='darkblue') + scale_y_continuous(labels=scales::percent) + xlab('# of cells') + ylab('% of molecules') + labs(subtitle='') + coord_cartesian(ylim=c(0,0.01))
	return(p)
}

#' @export
sc5mc.ontar_reads <- function(smat, intervals=NULL, return_orig=TRUE){
	if (!has_name(smat, 'tidy_cpgs')){
		loginfo("Getting tidy cpgs")
		smat <- smat.get_tidy_cpgs(smat, unique=TRUE)
	}
	
	if (!is.null(intervals)){
		reads <- smat$tidy_cpgs %>% distinct(read_id, num, .keep_all=T) %>% filter(end != '-')
		tcpgs <- gpatterns:::.gpatterns.filter_ontar_reads(reads, intervals, return_orig=return_orig)		
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
	rpu_fn <- glue('{smat$tidy_cpgs_dir}/stats/rpu.csv')
	if (file.exists(rpu_fn) && is.null(intervals)){
		rpu <- fread(rpu_fn) %>% as_tibble()
	} else {
		rpu <- sc5mc.reads_per_umi(smat, intervals)		
	}
	
	rpu <- rpu %>% 
		mutate(reads = cut(reads, breaks=c(0,1,4,10,1e6), include.lowest=TRUE, labels=c('1', '2-4', '5-10', '>10'))) %>% 
		group_by(reads) %>% 
		summarise(n = sum(n, na.rm=TRUE)) %>% 
		mutate(p = n / sum(n))	
	p <- rpu %>% ggplot(aes(x=reads, y=p)) + geom_col(width=0.7, fill='midnightblue') + scale_y_continuous(label=scales::percent) +  xlab('Reads / UMI') + ylab("% of reads") +theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5)) + labs(subtitle=glue("number of reads: {comify(sum(rpu$n))}"))
	return(p)	
}

sc5mc.plot_capture_reads_per_umi <- function(smat, intervals){
	rpu_offtar_fn <- glue('{smat$tidy_cpgs_dir}/stats/rpu_offtar.csv')
	if (file.exists(rpu_offtar_fn)){
		rpu_offtar <- fread(rpu_offtar_fn) %>% as_tibble()
	} else {
		rpu_offtar <- sc5mc.reads_per_umi(smat, gintervals.diff(gintervals.all(), intervals))
	}

	rpu_ontar_fn <- glue('{smat$tidy_cpgs_dir}/stats/rpu_ontar.csv')
	if (file.exists(rpu_ontar_fn)){
		rpu_ontar <- fread(rpu_ontar_fn) %>% as_tibble()
	} else {
		rpu_ontar <- sc5mc.reads_per_umi(smat, intervals)		
	}

	rpu <- rpu_offtar %>% 
		mutate(type = 'off target') %>% 
		bind_rows(rpu_ontar %>% mutate(type = 'on target')) %>% 
		mutate(reads = cut(reads, breaks=c(0,1,4,10,1e6), include.lowest=TRUE, labels=c('1', '2-4', '5-10', '>10'))) %>% 
		group_by(reads, type) %>% 
		summarise(n = sum(n, na.rm=TRUE)) %>% 
		group_by(type) %>%
		mutate(p = n / sum(n))	
	p <- rpu %>% ggplot(aes(x=reads, y=p, fill=type)) + geom_col(width=0.7, position='dodge') + scale_y_continuous(label=scales::percent) +  xlab('Reads / UMI') + ylab("% of reads") + ggsci::scale_fill_lancet(name='') + theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5)) + theme(legend.position="bottom")

	return(p)
}

sc5mc.reads_per_umi_insert_length <- function(smat, intervals=NULL){
	tcpgs <- sc5mc.ontar_reads(smat, intervals)	
	rpu <- tcpgs %>% mutate(reads = cut(num, breaks=c(0,1,4,10,1e6), include.lowest=TRUE, labels=c('1', '2-4', '5-10', '>10'))) %>% filter(!is.na(reads))
	return(rpu)
}

#' @export
sc5mc.plot_insert_length_reads_per_umi <- function(smat, intervals=NULL){
	rpu_fn <- glue('{smat$tidy_cpgs_dir}/stats/rpu_insert_length.csv')
	if (file.exists(rpu_fn) && is.null(intervals)){
		rpu <- fread(rpu_fn) %>% as_tibble()
	} else {
		rpu <- sc5mc.reads_per_umi_insert_length(smat, intervals)		
	}	
 	p <- rpu %>% ggplot(aes(x=abs(insert_len), color=reads)) + geom_density(size=0.9) + xlab('Insert length') +ylab('Density') + ggsci::scale_color_aaas(name = 'Reads / UMI') + theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5), legend.position = c(0.8, 0.8), legend.box.background = element_rect()) + coord_cartesian(xlim=c(0,450)) + guides(color=guide_legend(keywidth=0.5, keyheight=0.5))
	
	return(p)	
}


#' @export
sc5mc.frac_on_target <- function(smat, intervals){
	if (!has_name(smat, 'tidy_cpgs')){
		loginfo("Getting tidy cpgs")
		smat <- smat.get_tidy_cpgs(smat, unique=TRUE)
	}

	reads <- smat$tidy_cpgs %>% distinct(read_id, .keep_all=T) %>% filter(end != '-')
	ontar <- .gpatterns.filter_ontar_reads(reads, intervals)	
	stats <- reads %>% group_by(cell_id) %>% summarise(all_reads = sum(num)) %>% left_join(ontar %>%  group_by(cell_id) %>% summarise(ontar_reads = sum(num)), by='cell_id') %>% replace_na(ontar_reads = 0) %>% mutate(frac_ontar = ontar_reads / all_reads)
	
	return(stats)	
}

#' @export
sc5mc.plot_frac_on_target <- function(smat, intervals, min_reads=1e3){
	ontar_stats_fn <- glue('{smat$tidy_cpgs_dir}/stats/ontar_stats.csv')
	if (file.exists(ontar_stats_fn)){
		stats <- fread(ontar_stats_fn)
	} else {
		stats <- sc5mc.frac_on_target(smat, intervals)	
	}
	
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
	reg_stats_fn <- glue('{smat$tidy_cpgs_dir}/stats/regions_stats.csv')
	if(file.exists(reg_stats_fn)){
		reg_stats <- fread(reg_stats_fn) %>% as_tibble()
	} else {
		reg_stats <- sc5mc.reads_per_region(smat, intervals)		
	}
	
	reg_stats_p <- reg_stats %>% 
		ungroup() %>% 
		count(cov) %>% 
		mutate(cum = rev(cumsum(rev(n))))	
	
	p <- reg_stats_p %>% 
		filter(cov >= 5) %>%
		ggplot(aes(x=cov, y=cum)) + 
			geom_point(size=0.5, shape=19) + 			
			xlab("Reads") + 
			ylab("Regions with at least x reads") + 
			scale_y_continuous(labels=comify) +
			labs(subtitle=glue("number of reads: {comify(sum(reg_stats$cov))}"))
	return(p) 
}


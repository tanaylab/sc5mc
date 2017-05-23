#' @export
sc5mc.qc_plot <- function(smat, min_cpgs=1, max_cpgs=Inf, min_cells=1, max_cells=Inf){
	if (any(min_cpgs != 1, max_cpgs != Inf,  min_cells != 1, max_cells != Inf)){
		smat <- smat.filter_by_cov(smat, min_cells=min_cells, max_cells=max_cells, max_cpgs=max_cpgs, min_cpgs=min_cpgs)
	}
	cgp_mars <- sc5mc.plot_cpg_marginals(smat, type='percent')
	cell_mars <- sc5mc.plot_cell_marginals(smat, type='percent')
	pairs_mars <- sc5mc.plot_pairs_coverage(smat)
	p <- cowplot::plot_grid(cgp_mars, cell_mars, pairs_mars, align='hv', labels=LETTERS[1:3])
	return(p)
}

#' @export
sc5mc.plot_cpg_marginals <- function(smat, type='percent'){
	mars <- smat.cpg_marginals(smat) %>% filter(cells > 0) %>% select(cells)
	mars <- mars %>% group_by(cells) %>% summarise(cpgs = n()) %>% arrange(-cells) %>% mutate(c_cpgs = cumsum(cpgs), p = c_cpgs / sum(cpgs))
	x_breaks <- c(1,2,3,4,10,as.numeric(round(quantile(mars$cells, probs=c(0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9, 1))))) %>% unique %>% sort
	if (type == 'percent'){
		p <- mars %>% ggplot(aes(x=cells, y=p)) + geom_point(size=0.3) + scale_x_log10(breaks=x_breaks) + scale_y_continuous(label=scales::percent) + ylab('% of CpGs with at least x cells') + xlab('# of cells')	
	} else {
		p <- mars %>% ggplot(aes(x=cells, y=c_cpgs)) + geom_point(size=0.3) + scale_x_log10(breaks=x_breaks) + scale_y_log10(label=comify) + ylab('# of CpGs with at least x cells') + xlab('# of cells')
	}
	p <- p + labs(subtitle = qq('number of cells = @{comify(ncol(smat$cov))}, number of cpgs = @{comify(nrow(smat$cov))}'))
	return(p)	
}

#' @export
sc5mc.plot_cell_marginals <- function(smat, type='percent'){
	mars <- smat.cell_marginals(smat) %>% filter(cov > 0) %>% select(cov)
	mars <- mars %>% group_by(cov) %>% summarise(cells = n())  %>% arrange(-cov) %>% mutate(c_cells = cumsum(cells), p = c_cells / sum(cells))
	
	if (type == 'percent'){
		p <-mars %>% ggplot(aes(x=cov, y=p)) + geom_point(size=0.3) + scale_x_log10(label=comify) + scale_y_continuous(label=scales::percent) + ylab('% of cells with at least x CpGs') + xlab('# of CpGs')
	} else {
		p <- mars %>% ggplot(aes(x=cov, y=c_cells)) + geom_point(size=0.3) + scale_x_log10(label=comify) + scale_y_continuous(label=comify) + ylab('# of cells with at least x CpGs') + xlab('# of CpGs')
	}

	p <- p + labs(subtitle = qq('number of cells = @{comify(ncol(smat$cov))}, number of cpgs = @{comify(nrow(smat$cov))}'))
	return(p)	
}

#' @export
sc5mc.plot_pairs_coverage <- function(smat){	
	mars <- smat.cell_pairs_marginals(smat)
	cmars <- mars %>% count(ntot) %>% arrange(-ntot) %>% mutate(n = n / 2, cn = cumsum(n), p = cn / sum(n))

	x_breaks <- as.numeric(round(quantile(cmars$ntot, probs=c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1))))
	p <- cmars %>% ggplot(aes(x=ntot, y=p)) + geom_point(size=0.3) + xlab('# of common CpGs') + ylab('% of cell pairs') + scale_x_log10(label=comify, breaks=x_breaks) + scale_y_continuous(labels=scales::percent) + labs(subtitle = qq('number of cells = @{comify(length(unique(mars$cell1)))}, number of pairs = @{comify(nrow(mars) / 2)}'))
	return(p)
	
}
#' @import Matrix
#' @import gpatterns


tidy2smat <- function(data, row, column, value, ...){
    row_u <- unique(data[[row]])
    i <- match(data[[row]], row_u)

    col_u <- unique(data[[column]])
    j <- match(data[[column]], col_u)

    val <- data[[value]]
    
    Matrix::sparseMatrix(i = i, j = j, x = val, dimnames = list(row_u, col_u), ...)
}

.tidy_calls2smat <- function(tidy_calls){
    smat <- alply(c('meth', 'unmeth', 'cov'), 1, function(x) {
        message(sprintf('creating %s', x))
        tidy2smat(tidy_calls, 'coord', 'track', x)                    
    }, .parallel=TRUE)

    names(smat) <- c('meth', 'unmeth', 'cov')
  
    message('creating intervs')
    smat$intervs <- tibble(coord = rownames(smat[[1]])) %>% separate(coord, c('chrom', 'start', 'end'))

    return(smat)
}

tcpgs2calls <- function(tcpgs, track){
   tcpgs %>% 
        arrange(chrom, cg_pos) %>% 
        distinct(chrom, cg_pos, .keep_all=T) %>% 
        mutate(start=cg_pos, end=start+1, track=track) %>% 
        select(track, chrom, start, end, meth) %>% 
        unite('coord', chrom:end)
}

#' @export
smat.import_from_bams <- function(metadata, groot, prefix=NULL, workdir=tempdir(), ...){
    bam2calls <- function(bams, lib, workdir=workdir, ...){
        track_workdir <- tempfile(tmpdir=workdir)
        system(sprintf('mkdir -p %s', track_workdir))
        on.exit(system(sprintf('rm -rf %s', track_workdir)))
        gpatterns.import_from_bam(bams, workdir=track_workdir, steps=c('bam2tidy_cpgs', 'filter_dups'), groot=groot, ...)

        tcpgs_dir <- paste0(track_workdir, '/tidy_cpgs_uniq')
        tcpgs <- .gpatterns.get_tidy_cpgs_from_dir(tcpgs_dir)
        
        if (is.null(tcpgs)){
            return(NULL)
        }
        calls <- tcpgs2calls(tcpgs, lib)

        return(calls)
    }

    gsetroot(groot)
    
    commands <- plyr::daply(metadata, .(lib), function(x) {
        bams <- paste(qqv("'@{x$bam}'"), collapse=', ')
        qq("bam2calls(c(@{bams}), lib = '@{x$lib[1]}', workdir='@{workdir}', ...)")
    })

    res <- gcluster.run2(command_list=commands, io_saturation=T)

    tidy_calls <- res %>% map('retv') %>% compact() %>% map_df(~ .x) 
    tidy_calls <- tidy_calls %>% mutate(unmeth = if_else(meth == 0, 1, 0), cov=1)

    smat <- .tidy_calls2smat(tidy_calls)

    if (!is.null(prefix)){
        smat.save(smat, prefix)
    }
    return(smat)    
}

#' @export
smat.import_from_tracks <- function(tracks, libs, prefix=NULL, threads=5){
    get_tidy_meth_calls <- function(track, uniq=TRUE, dsn=NULL) {        
        tcpgs <- gpatterns::gpatterns.get_tidy_cpgs(track)    
        if (is.null(tcpgs)){
            return(NULL)
        }    
        tcpgs2calls(tcpgs, track)     
    }
    commands <- paste0('get_tidy_meth_calls(tracks[', 1:length(tracks), '])')

    res <- gcluster.run2(command_list=commands, threads=threads)
    
    tidy_calls <- res %>% map('retv') %>% compact() %>% map_df(~ .x) 
    tidy_calls <- tidy_calls %>% mutate(unmeth = if_else(meth == 0, 1, 0), cov=1)

    smat <- .tidy_calls2smat(tidy_calls)

    if (!is.null(prefix)){
        smat.save(smat, prefix)
    }
    return(smat)
}

#' @export
smat.save <- function(smat, prefix){
    res <- alply(c('meth', 'unmeth', 'cov'), 1, function(x) {
        message(sprintf('creating %s', x))        
        Matrix::writeMM(smat[[x]], paste0(prefix, '.', x))        
    }, .parallel=TRUE)

    tibble(track = colnames(smat$cov)) %>% fwrite(paste0(prefix, '_colnames.tsv'))

    message('creating intervs')
    tibble(coord = rownames(smat$cov)) %>% separate(coord, c('chrom', 'start', 'end')) %>% fwrite(paste0(prefix, '_intervs.tsv'))
    invisible(smat)   
}

#' @export
smat.load <- function(conf){      
    mat_colnames <- read_csv(conf$sparse_matrix$colnames, col_types=cols(track = col_character())) %>% .$track

    smat <- alply(c('cov', 'meth', 'unmeth'), 1, function(x) {
    	message(sprintf('doing %s', x))
    	m <- readMM(conf$sparse_matrix[[x]]) * 1
    	colnames(m) <- mat_colnames
    	m}, .parallel=TRUE)

    names(smat) <- c('cov', 'meth', 'unmeth')

    obs_cpgs <- fread(conf$sparse_matrix$intervs) %>% mutate(id = 1:n()) %>% tbl_df
    smat$intervs <- obs_cpgs
    return(smat)
}

#' @export
smat.cpg_marginals <- function(smat){
	mars <- rowSums(smat$cov)
	return(smat$intervs %>% select(chrom, start, end) %>% mutate(cells = mars))
}

#' @export
smat.cell_marginals <- function(smat){
	mars <- colSums(smat$cov)
	return(tibble(track=colnames(smat$cov), cov=mars))
}

#' @export
smat.cell_pairs_marginals <- function(smat){
    ntot <- crossprod(smat$cov) %>% as.matrix() %>% reshape2::melt() %>% rename(cell1=Var1, cell2=Var2, ntot=value)
    return(ntot)
}

#' @export
smat.filter <- function(smat, intervs=NULL, cols=NULL, ...){
	if (!is.null(intervs)){
		intervs <- smat$intervs %>% 
    	gintervals.intersect(intervs)     		
	}     
    
    return(smat.filter_cpgs(smat, intervs, cols=cols, ...))    
}

#' @export
smat.filter_cpgs <- function(smat, intervs=NULL, cols=NULL){
	if (is.null(intervs)){
		new_mat_intervs <- smat$intervs		
		ids <- 1:nrow(smat$cov) 	
	} else {
		new_mat_intervs <- smat$intervs %>% inner_join(intervs, by=c('chrom', 'start', 'end')) %>% arrange(id)		
		ids <- new_mat_intervs$id 	
	}
	
		
	if (is.null(cols)){
		cols <- colnames(smat$cov)
	}	

	new_smat <- smat
	new_smat$cov <- smat$cov[ids, cols]
    new_smat$meth <- smat$meth[ids, cols]
    new_smat$unmeth <- smat$unmeth[ids, cols]
    new_smat$intervs <- new_mat_intervs %>% mutate(id = 1:n())
    return(new_smat)
}

#' @export
smat.filter_by_cov <- function(smat, min_cpgs=1, max_cpgs=Inf, min_cells=1, max_cells=Inf){
	cell_filter <- which(between(colSums(smat$cov), min_cpgs, max_cpgs))
	cg_filter <- which(between(rowSums(smat$cov), min_cells, max_cells))

	new_smat <- smat
	new_smat$cov <- smat$cov[cg_filter, cell_filter]
    new_smat$meth <- smat$meth[cg_filter, cell_filter]
    new_smat$unmeth <- smat$unmeth[cg_filter, cell_filter]
    new_smat$intervs <- smat$intervs[cg_filter, ] %>% mutate(id = 1:n())

    return(new_smat)
}


smat.summarise_by_track <- function(smat, track, breaks, include.lowest=TRUE, group_name='group'){
    opt <- getOption('gmax.data.size')
    options(gmax.data.size=1e9)
    on.exit(options(gmax.data.size=opt))
    message(qq('extracting @{track}'))
    track_df <- gextract(track, intervals=smat$intervs, iterator=smat$intervs, colnames='group') %>% arrange(intervalID)
    groups <- cut(track_df$group, breaks=breaks, include.lowest=include.lowest) 
    message(qq('summarising...'))
    res <- smat.summarise(smat, groups=groups, group_name=group_name)
    return(res)
}

#' @export
smat.summarise <- function(smat, groups, group_name='group'){
    res <- smat$intervs %>% 
        mutate(group = groups) %>% 
        ddply(.(group), function(x) broom::tidy(smat$meth[x$id, ]) %>% 
            group_by(column) %>% 
            summarise(ncpgs = n(), meth=mean(value)),
             .parallel=T) %>% 
        tbl_df %>% 
        set_names(c(group_name, 'track', 'ncpgs', 'meth'))
    
    return(res)
}

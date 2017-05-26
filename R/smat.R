#' @import Matrix
#' @import gpatterns

## Import functions
tidy2smat <- function(data, row, column, value, ...){
    row_u <- unique(data[[row]])
    i <- match(data[[row]], row_u)

    col_u <- unique(data[[column]])
    j <- match(data[[column]], col_u)

    val <- data[[value]]

    Matrix::sparseMatrix(i = i, j = j, x = val, dimnames = list(row_u, col_u), ...)
}

.tidy_calls2smat <- function(tidy_calls){
    smat <- plyr::alply(c('meth', 'unmeth', 'cov'), 1, function(x) {
        message(sprintf('creating %s', x))
        tidy2smat(tidy_calls, 'coord', 'track', x)
    }, .parallel=TRUE)

    names(smat) <- c('meth', 'unmeth', 'cov')

    message('creating intervs')
    smat$intervs <- tibble(coord = rownames(smat[[1]])) %>% 
        separate(coord, c('chrom', 'start', 'end')) %>% mutate(id = 1:n())

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

#' Create smat object from mapped reads (bam files)
#'
#' @param metadata data frame with 'lib' column with the library name (ususally cell id)
#' and 'bam' column with bam files. Multiple bams of the same cell should appear
#' in separate rows.
#'
#' @param groot misha db path
#' @param prefix path prefix in which to save the smat files. if NULL the object
#' would not be written to disk
#' @param workdir temporary directory to use for bam parsing
#' @param use_sge use sun grid engine cluster
#' @param ... additional parameters to \code{gpatterns::gpatterns.import_from_bam}
#'
#' @export
smat.from_bams <- function(metadata, groot, prefix=NULL, workdir=tempdir(), use_sge = TRUE, ...){
    bam2calls <- function(bams, lib, workdir=workdir, use_sge=TRUE, ...){
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
    if (use_sge){
        res <- gcluster.run2(command_list=commands, io_saturation=T)
        tidy_calls <- res %>% map('retv') %>% compact() %>% map_df(~ .x)
    } else {
        tidy_calls <- map(commands, ~ eval(parse(text=.x))) %>% compact() %>% map_df(~ .x)
    }

    tidy_calls <- tidy_calls %>% mutate(unmeth = if_else(meth == 0, 1, 0), cov=1)

    smat <- .tidy_calls2smat(tidy_calls)

    if (!is.null(prefix)){
        smat.save(smat, prefix)
    }
    return(smat)
}

#' Create smat object from misha tracks
#'
#' @param tracks names of tracks
#'
#' @param libs names of new object smat columns
#' @param prefix path prefix in which to save the smat files. if NULL the object
#' would not be written to disk
#' @param use_sge use sun grid engine cluster
#' @param ... additional parameters to gcluster.run2
#'
#' @return smat object
#'
#' @export
smat.from_tracks <- function(tracks, libs, prefix=NULL, use_sge=TRUE, ...){
    get_tidy_meth_calls <- function(track, uniq=TRUE, dsn=NULL) {
        tcpgs <- gpatterns::gpatterns.get_tidy_cpgs(track)
        if (is.null(tcpgs)){
            return(NULL)
        }
        tcpgs2calls(tcpgs, track)
    }
    commands <- paste0('get_tidy_meth_calls(tracks[', 1:length(tracks), '])')
    if (use_sge){
        res <- gcluster.run2(command_list=commands, ...)
        tidy_calls <- res %>% map('retv') %>% compact() %>% map_df(~ .x)
    } else {
        tidy_calls <- map(commands, ~ eval(parse(text=.x))) %>% compact() %>% map_df(~ .x)
    }


    tidy_calls <- tidy_calls %>% mutate(unmeth = if_else(meth == 0, 1, 0), cov=1)

    smat <- .tidy_calls2smat(tidy_calls)

    if (!is.null(prefix)){
        smat.save(smat, prefix)
    }
    return(smat)
}

#' Create smat object from data frame
#' @param df intervals (chrom, start, end fields) data frame with additonal fields:
#' meth (methylated calls), unmeth (unmethylated calls) and cov (total coverage)
#'
#' @return smat object
#' @export
smat.from_df <- function(df){
    .tidy_calls2smat(df %>% unite('coord', chrom, start, end))
}

#' export smat object to data frame
#'
#' @param smat smat object
#'
#' @param coords if TRUE - chrom,start,end fields are united to a single field named 'coords'
#'
#' @return intervals (chrom, start, end fields) data frame with additonal fields:
#' meth (methylated calls), unmeth (unmethylated calls) and cov (total coverage)
#'
#' @export
smat.to_df <- function(smat, coords=TRUE){
    tmeth <- broom::tidy(smat$meth) %>% tbl_df %>% rename(id=row, cell=column, meth=value)
    tunmeth <- broom::tidy(smat$unmeth) %>% tbl_df %>% rename(id=row, cell=column, unmeth=value)
    tidy <- tmeth %>% mutate(unmeth = tunmeth[['unmeth']], cov = meth + unmeth)
    if (coords){
        tidy <- tidy %>% separate(id, c('chrom', 'start', 'end'))
    }
    return(tidy)
}

#' Save smat object to disk
#'
#' @param smat smat object
#'
#' @param prefix path prefix of the location of the smat files
#'
#' @return inivisibly returns the smat object
#'
#' @export
smat.save <- function(smat, prefix){
    res <- plyr::alply(c('meth', 'unmeth', 'cov'), 1, function(x) {
        message(sprintf('creating %s', x))
        Matrix::writeMM(smat[[x]], paste0(prefix, '.', x))
    }, .parallel=TRUE)

    tibble(cell = smat.colnames(smat)) %>% fwrite(paste0(prefix, '_colnames.tsv'))

    message('creating intervs')
    tibble(coord = rownames(smat$cov)) %>% separate(coord, c('chrom', 'start', 'end')) %>% fwrite(paste0(prefix, '_intervs.tsv'))
    invisible(smat)
}


#' Load smat object from disk
#'
#' @param prefix path prefix of the location of the smat files
#'
#' @return smat object
#'
#' @export
smat.load <- function(prefix){
    conf <- list(sparse_matrix=list(
        meth = paste0(prefix, '.meth'),
        unmeth = paste0(prefix, '.unmeth'),
        cov = paste0(prefix, '.cov'),
        intervs = paste0(prefix, '_intervs.tsv'),
        colnames = paste0(prefix, '_colnames.tsv')
        ))
    return(smat.from_conf(conf))
}

#' Load smat object from disk
#'
#' @param conf list with 'sparse_matrices' field that holds a list with:
#' cov - path of coverage sparse matrix
#' meth - path of meth sparse matrix
#' unmeth - path of unmeth sparse matrix
#' intervs - path of intervs file
#' colnames - path of colnames file
#'
#' @return smat object
#'
#' @export
smat.from_conf <- function(conf){
    mat_colnames <- read_csv(conf$sparse_matrix$colnames, col_types=cols(col_character()))[[1]]

    smat <- plyr::alply(c('cov', 'meth', 'unmeth'), 1, function(x) {
        message(sprintf('doing %s', x))
        m <- readMM(conf$sparse_matrix[[x]]) * 1
        colnames(m) <- mat_colnames
        return(m)
          }, .parallel=TRUE)

    names(smat) <- c('cov', 'meth', 'unmeth')

    obs_cpgs <- fread(conf$sparse_matrix$intervs) %>% mutate(id = 1:n()) %>% tbl_df
    smat$intervs <- obs_cpgs
    return(smat)
}

#' Print smat object
#' 
#' @param smat smat object
#' 
#' @export
smat.info <- function(smat){
    ncells <- comify(ncol(smat[['cov']]))
    ncpgs <- comify(nrow(smat[['cov']]))
    attrs <- paste(names(smat), collapse=', ')
    cat(qq('smat object\n@{ncpgs} CpGs X @{ncells} cells\n'))
    cat(qq('attributes: @{attrs}\n'))
}

#' Get colnames of smat object
#' 
#' @param smat smat object
#' 
#' @export
smat.colnames <- function(smat){
    return(colnames(smat[['cov']]))
}


#' convert column names to ids
cols2ids <- function(smat, cols=NULL){
    if (!is.null(cols)){
        ids <- which(smat.colnames(smat) %in% cols)
    } else {
        ids <- 1:ncol(smat$cov)
    }
    return(ids)
}


#' calculate marginal coverage of CpGs
#'
#' @param smat smat object
#' @param cols column names of columns to use. If NULL all columns (cells) would be returned
#'
#' @return intervals set (chrom, start, end) with the cpgs and their marginal
#' coverage (cells field)
#' 
#' @export
smat.cpg_marginals <- function(smat, cols=NULL){	
    ids <- cols2ids(smat, cols)

    mars <- rowSums(smat$cov[, ids])

	return(smat$intervs %>% select(chrom, start, end) %>% mutate(cells = mars))
}

#' calculate marginal coverage of cells
#'
#' @param smat smat object
#' @param cols column names of columns to use. If NULL all columns (cells) would be returned
#'
#' @return tibble with cells ('cell' field) and their marignal coverage ('cov' field)
#' @export
smat.cell_marginals <- function(smat, cols=NULL){
    ids <- cols2ids(smat, cols)
	mars <- colSums(smat$cov[, ids])

    return(tibble(cell=smat.colnames(smat)[ids], cov=mars))
}

#' Calculate the joint coverage of all pairs of cells
#'
#' @param smat smat object
#' @param cols column names of columns to use. If NULL all columns (cells) would be returned
#' 
#' @return tibble with cell pairs (cell1, cell2 fields) and 'ntot' with the
#' number of jointly covered CpGs
#'
#' @export
smat.cell_pairs_marginals <- function(smat, cols=NULL){
    ids <- cols2ids(smat, cols)   

    ntot <- crossprod(smat$cov[, ids]) %>% as.matrix() %>% reshape2::melt() %>% rename(cell1=Var1, cell2=Var2, ntot=value)

    ntot <- combn(smat.colnames(smat)[ids], 2) %>% t() %>% tbl_df %>% set_names(c('cell1', 'cell2')) %>% inner_join(ntot, by=c('cell1', 'cell2'))
    
    return(ntot)
}

#' Calculate the joint coverage of all pairs of CpGs
#'
#' @param smat smat object
#' 
#' @param min_cells minimal number of cells per CpG
#' @param cols column names of columns to use. If NULL all columns (cells) would be returned
#' 
#' @return tibble with CpG coverage ('ncells') and the number of CpG pairs with that number 
#' jointly covered cells ('npairs')
#'
#' @export
smat.cpg_pairs_marginals <- function(smat, min_cells=30, cols=NULL){
    ids <- cols2ids(smat, cols)

    cgs <- which(rowSums(smat$cov[, ids]) >= min_cells)
    
    ntot <- tcrossprod(smat$cov[cgs, ids]) %>% broom::tidy()
    # add explicit zeroes and remove double counting
    res <-  combn(1:length(cgs), 2) %>% t() %>% tbl_df %>% set_names(c('row', 'column')) %>% left_join(ntot, by=c('row', 'column')) %>% mutate(value = if_else(is.na(value), 0, value))
    res <- res %>% count(value) %>% rename(ncells=value, npairs=n)

    return(res)
}

#' Filter smat object
#'
#' @param smat smat object
#'
#' @param intervs intervals set (df with chrom,start,end fields) to filter smat by.
#' only CpGs covered by the intervals would be returned. if NULL all CpGs would be returned.
#' @param cols column names of columns to keep. If NULL all columns (cells) would be returned
#' @param ids vector ids of CpGs to filter by.
#'
#' @return smat object with CpGs covered by intervs, and cells that in cols / ids.
#'
#' @export
smat.filter <- function(smat, intervs=NULL, cols=NULL, ids=NULL){
	if (!is.null(intervs)){
		intervs <- smat$intervs %>%
    	gintervals.intersect(intervs)
	}

    return(smat.filter_cpgs(smat, cpg_intervs=intervs, cols=cols, ids=ids))
}

#' Select cells from smat object
#' 
#' @inheritParams smat.filter
#' @export
smat.select <- function(smat, cols=NULL){
    smat.filter(smat, cols=cols)
}

#' Filter smat object by CpG intervals
#'
#' @param smat smat object
#'
#' @param cpg_intervs intervals set of CpGs (df with chrom,start,end fields) to filter smat by.
#' only these CpGs would be returned. if NULL all CpGs would be returned.
#' @param cols column names of columns to keep. If NULL all columns (cells) would be returned
#' @param ids vector ids of CpGs to filter by.
#'
#' @return smat object with CpGs covered by intervs, and cells that in cols / ids.
#'
#' @export
smat.filter_cpgs <- function(smat, cpg_intervs=NULL, cols=NULL, ids=NULL){
    if (!is.null(ids)){
        new_mat_intervs <- smat$cpg_intervs %>% filter(id %in% ids)
    } else if (!is.null(cpg_intervs)) {
        new_mat_intervs <- smat$cpg_intervs %>% 
            inner_join(cpg_intervs, by=c('chrom', 'start', 'end')) %>% 
            arrange(id)
        ids <- new_mat_intervs$id
    } else {
        new_mat_intervs <- smat$cpg_intervs
        ids <- 1:nrow(smat$cov)
    }

	if (is.null(cols)){
		cols <- smat.colnames(smat)
	} else {
        cols <- which(smat.colnames(smat) %in% cols)
    }

	new_smat <- smat
	new_smat$cov <- smat$cov[ids, cols]
    new_smat$meth <- smat$meth[ids, cols]
    new_smat$unmeth <- smat$unmeth[ids, cols]
    new_smat$intervs <- new_mat_intervs %>% mutate(id = 1:n())
    return(new_smat)
}

#' Filter smat object by cpgs / cells coverage
#'
#' @param smat smat object
#'
#' @param min_cpgs minimal number of CpGs
#' @param max_cpgs maximal number of CpGs
#' @param min_cells minimal number of cells
#' @param max_cells maximal number of cells
#'
#' @export
smat.filter_by_cov <- function(smat, min_cpgs=1, max_cpgs=Inf, min_cells=1, max_cells=Inf){
    if (all(min_cpgs == 1, max_cpgs == Inf, min_cells == 1, max_cells == Inf)){
        return(smat)
    }
	cell_filter <- which(between(colSums(smat$cov), min_cpgs, max_cpgs))
	cg_filter <- which(between(rowSums(smat$cov), min_cells, max_cells))

    if (length(cg_filter) == 0){
        stop('no cpgs match the criteria')
    }
    if (length(cell_filter) == 0){
        stop('no cells match the criteria')
    }

	new_smat <- smat
	new_smat$cov <- smat$cov[cg_filter, cell_filter]
    new_smat$meth <- smat$meth[cg_filter, cell_filter]
    new_smat$unmeth <- smat$unmeth[cg_filter, cell_filter]
    new_smat$intervs <- smat$intervs[cg_filter, ] %>% mutate(id = 1:n())

    return(new_smat)
}

#' Group CpGs by a misha track, and summarise per cell
#'
#' @param smat smat object
#' @param track misha track expression to summarise by
#' @param breaks breaks that determine the bins of the track expression
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param group_name name of the column with the group in the returned data frame
#' @param parallel compute parallely per group (using doMC package)
#'
#' @return data frame with the following fields:
#' 'group_name', 'cell', 'ncpgs', 'meth', 'unmeth', 'avg'
#'
#' @export
smat.summarise_by_track <- function(smat, track, breaks, include.lowest=TRUE, group_name='group', parallel=TRUE){
    opt <- getOption('gmax.data.size')
    options(gmax.data.size=1e9)
    on.exit(options(gmax.data.size=opt))
    message(qq('extracting @{track}'))
    track_df <- gextract(track, intervals=smat$intervs, iterator=smat$intervs, colnames='group') %>% arrange(intervalID)
    groups <- cut(track_df$group, breaks=breaks, include.lowest=include.lowest)
    message(qq('summarising per group...'))
    res <- smat.summarise(smat, groups=groups, group_name=group_name, parallel=parallel)
    return(res)
}

#' Group CpGs by intervals, and summarise per cell
#'
#' @param smat smat object
#'
#' @param intervals intervals set. CpGs covered by each interval would be summarised
#' @param return_smat return smat object instead of data frame
#'
#' @return  if \code{return_smat} - smat object with rows as the intervals.
#' if not - data frame with the following fields:
#' 'chrom', 'start', 'end', 'cell', 'ncpgs', 'cov', 'meth', 'unmeth', 'avg'
#'
#' @export
smat.summarise_by_intervals <- function(smat, intervals, return_smat=FALSE){
    opt <- getOption('gmax.data.size')
    options(gmax.data.size=1e9)
    on.exit(options(gmax.data.size=opt))

    if ('strand' %in% colnames(intervals)){
        intervals1 <- intervals %>% select(-strand)
    } else {
        intervals1 <- intervals
    }
    message('finding overlapping CpGs')
    neighbours <- smat$intervs %>% gintervals.neighbors1(intervals1, na.if.notfound=FALSE, mindist=0, maxdist=0)

    smat_f <- smat.filter(smat, ids=neighbours$id)

    message(qq('summarising per group...'))
    tcpgs <- smat.to_df(smat_f, coords=FALSE)

    res <- tcpgs %>%
        inner_join(neighbours %>% select(id, chrom=chrom1, start=start1, end=end1), by='id') %>%
        group_by(chrom, start, end, cell) %>%
        summarise(ncpgs=n(), cov = sum(cov), meth = sum(meth), unmeth = sum(unmeth), avg = meth / cov) %>%
        ungroup

    if (return_smat){
        res <- smat.from_df(res)
    }

    return(res)
}

#' Summarise groups of CpGs per cell
#'
#' @param smat sc5mc smat object
#' @param groups character / factor vector with a group per CpG
#' @param group_name column name of the group
#' @param parallel compute parallely per group (using doMC package)
#'
#'
#' @return data frame with the following fields:
#' 'group_name', 'cell', 'ncpgs', 'meth', 'unmeth', 'avg'
#' @export
smat.summarise <- function(smat, groups, group_name='group', parallel=TRUE){
    summarise_per_group <- function(x){
        if (nrow(x) == 1){
            m <- smat$meth[x$id, ]
            um <- smat$unmeth[x$id, ]
            cov <- smat$cov[x$id, ]
            avg <- m / cov
        } else {
            m <- colSums(smat$meth[x$id, ])
            um <- colSums(smat$unmeth[x$id, ])
            cov <- colSums(smat$cov[x$id, ])
            avg <- m / cov
        }

        tibble(cell=names(cov), ncpgs = cov, meth=m, unmeth=um, avg=avg)
    }

    res <- smat$intervs %>%
        mutate(group = groups) %>%
        plyr::ddply(.(group), summarise_per_group, .parallel=parallel) %>%
        tbl_df %>%
        filter(ncpgs > 0) %>%
        set_names(c(group_name, 'cell', 'ncpgs', 'meth', 'unmeth', 'avg'))

    return(res)
}

#' Summarises CpGs by groups pf cells
#'
#' @param smat sc5mc smat object
#' @param groups_df data frame with column 'cell' with cell names, and additional columns with groups. if NULL CpGs from all cells would be summarised.
#' @param groups columns of group_df to group by.
#' @param min_cells mininmal number of cells per CpG
#' @param max_cells maximal number of cells cells per CpG
#'
#' @return tidy data frame with summariy statistics for each group and CpG
#'
#' @export
smat.summarise_cpgs <- function(smat, groups_df=NULL, groups=NULL, min_cells=1, max_cells=Inf){
    summarise_per_group <- function(x){
        cells <- match(x$cell, smat.colnames(smat))
        cells <- cells[!is.na(cells)]
        cgs <- which(between(rowSums(smat[['cov']][, cells]), min_cells, max_cells))

        tibble(id = cgs,
               m = rowSums(smat[['meth']][cgs, cells]),
               n = rowSums(smat[['cov']][cgs, cells]))
    }

    rm_group <- FALSE
    if (is.null(groups_df) || is.null(groups)){
        groups_df <- tibble(cell = smat.colnames(smat)) %>% mutate(group = 'group')
        groups <- 'group'
        rm_group <- TRUE
    }
    res <- groups_df %>% select(one_of('cell', groups)) %>% group_by_(groups) %>% do(summarise_per_group(.)) %>% ungroup %>% mutate(unmeth = n - m, avg = m / n)
    res <- res %>% left_join(smat$intervs, by='id') %>% select(chrom, start, end, id, one_of(groups), cov=n, meth = m, unmeth, avg)

    if (rm_group){
        res <- res %>% select(-group)
    }

    return(res)
}

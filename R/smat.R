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
smat.to_df <- function(smat, coords=FALSE){    
    tmeth <- broom::tidy(smat$meth) %>% tbl_df %>% rename(id=row, cell=column, meth=value)
    tunmeth <- broom::tidy(smat$unmeth) %>% tbl_df %>% rename(id=row, cell=column, unmeth=value)
    tidy <- tmeth %>% mutate(unmeth = tunmeth[['unmeth']], cov = meth + unmeth)
    tidy <- tidy %>% mutate(id = as.character(id)) %>% left_join(smat$intervs %>% mutate(id = as.character(id)), by='id') %>% select(chrom, start, end, id, cell, meth, unmeth, cov)
    if (coords){
        tidy <- tidy %>% unite('coords', chrom, start, end)
    }
    return(tidy)
}

#' Create track from smat marginals
#' @export
smat.to_marginal_track <- function(smat, track, description, cols=NULL, overwrite=FALSE){
    pileup <-  smat.cpg_avg_marginals(smat, cols=cols) %>% select(chrom, start, end, cov, meth, unmeth)
    gpatterns:::.gpatterns.import_intervs_table(track, description, pileup, overwrite=overwrite)
}

#' join two smat objects
#' 
#' @param smat_r first smat object
#' @param smat_l seconf smat object
#' @param type join type. can be 'inner', 'full', 'left', 'right'. see xx_join by dplyr package.
#' @param prefix_l prefix for cell names from first smat in case of conflicts
#' @param prefix_r prefix for cell names from second smat in case of conflicts
#' 
#' @return smat object with intervals joined by \code{type}_join
#' 
#' @export
smat.join <- function(smat_r, smat_l, type='full', prefix_l='l_', prefix_r='r_'){ 
    join_func <- switch(type, 
        inner = gintervals.intersect, 
        full = gintervals.union, 
        left = function(l, r) gintervals.intersect(gintervals.union(l,r), l), 
        right = function(l, r) gintervals.intersect(gintervals.union(l,r), r))
    
    message(sprintf('merging intervals, join type: %s', type))
    intervs <- join_func(smat_l$intervs, smat_r$intervs) %>% mutate(id = 1:n())

    message('merging sparse matrices')
    df_r <- smat.to_df(smat_r) %>% select(-id)        
    df_l <- smat.to_df(smat_l) %>% select(-id)       
    
    conf_names <- df_l %>%
        distinct(cell) %>% 
        inner_join(df_r %>% 
            distinct(cell), by='cell') %>%
        pull(cell)

    if (length(conf_names) > 0){
        warning(qq('conflicting names, adding prefixes: @{paste(conf_names, collapse=", ")}'))
        df_l <- df_l %>% mutate(cell = paste0(prefix_l, cell))
        df_r <- df_r %>% mutate(cell = paste0(prefix_r, cell))
    }
    
    # smat.from_df
    df <- bind_rows(df_l, df_r) %>% inner_join(intervs, by=c('chrom', 'start', 'end'))    

    message('creating smat from df')
    smat <- smat.from_df(df, intervs=intervs)

    if (has_stats(smat_l) && has_stats(smat_r) && length(conf_names) == 0){
        smat$stats <- bind_rows(smat_l$stats, smat_r$stats)
    }
    if (has_cell_metadata(smat_l) && has_cell_metadata(smat_r) && length(conf_names) == 0){
        smat$cell_metadata <- bind_rows(smat_l$cell_metadata, smat_r$cell_metadata)
    }    

    gc()
    
    return(smat)    
}

#' join multiple smat objects
#' 
#' @param ... smat objects to join
#' @inheritParams smat.join
#' 
#' @return joined smat objects
#' 
#' @export
smat.multi_join <- function(..., type='full', prefix_l='l_', prefix_r='r_'){
    matrices <- list(...)    
    purrr::reduce(matrices, smat.join, type=type, prefix_l=prefix_l, prefix_r=prefix_r)
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
        message(sprintf('saving %s', x))
        Matrix::writeMM(smat[[x]], paste0(prefix, '.', x))
    }, .parallel=TRUE)

    tibble(cell = smat.colnames(smat)) %>% fwrite(paste0(prefix, '_colnames.tsv'), sep='\t')

    message('saving intervs')
    fwrite(smat$intervs, paste0(prefix, '_intervs.tsv'), sep='\t')    

    if (has_stats(smat)){
        message('saving stats')
        fwrite(smat$stats, paste0(prefix, '_stats.tsv'), sep='\t')
    }

    if (has_name(smat, 'cell_metadata')){
        message('saving cell metadata')
        fwrite(smat$cell_metadata, paste0(prefix, '_cell_metadata.tsv'), sep='\t')
    }

    attributes <- list()
    for (attr in c('name', 'description')){
        if (has_name(smat, attr)){
            attributes[[attr]] <- smat[[attr]]    
        }        
    }

    if (length(attributes) > 0){
        message('saving attributes')
        readr::write_lines(yaml::as.yaml(attributes), paste0(prefix, '_attributes.yaml'))
    }

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

    if (file.exists(paste0(prefix, '_stats.tsv'))){
        conf$sparse_matrix$stats <- paste0(prefix, '_stats.tsv')
    }

    attributes_file <- paste0(prefix, '_attributes.yaml') 
    if (file.exists(attributes_file)){
        conf$sparse_matrix$attributes_file <- attributes_file        
    }

    cell_metadata_file <- paste0(prefix, '_cell_metadata.tsv')
    if (file.exists(cell_metadata_file)){
        conf$sparse_matrix$cell_metadata <- cell_metadata_file
    }
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
#' attributes_file - path of attributes yaml file (optional)
#'
#' @return smat object
#'
#' @export
smat.from_conf <- function(conf){
    mat_colnames <- fread(conf$sparse_matrix$colnames) %>% as.tibble() %>% pull(1)

    smat <- plyr::alply(c('cov', 'meth', 'unmeth'), 1, function(x) {
        message(sprintf('loading %s', x))
        m <- readMM(conf$sparse_matrix[[x]]) * 1
        colnames(m) <- mat_colnames
        return(m)
          }, .parallel=TRUE)

    names(smat) <- c('cov', 'meth', 'unmeth')

    obs_cpgs <- fread(conf$sparse_matrix$intervs) %>% mutate(id = 1:n()) %>% tbl_df()
    smat$intervs <- obs_cpgs

    if (has_name(conf$sparse_matrix, 'stats')){
        smat$stats <- fread(conf$sparse_matrix$stats) %>% tbl_df()
    }
    if (has_name(conf$sparse_matrix, 'attributes_file')){
        attributes <- yaml::yaml.load_file(conf$sparse_matrix$attributes_file)
        for (attr in names(attributes)){
            smat[[attr]] <- attributes[[attr]]
        }
    }
    if (has_name(conf$sparse_matrix, 'cell_metadata')){
        smat$cell_metadata <- fread(conf$sparse_matrix$cell_metadata) %>% tbl_df()
    }
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

#' calculate marginal average methylation of CpGs
#'
#' @param smat smat object
#' @param cols column names of columns to use. If NULL all columns (cells) would be returned
#'
#' @return intervals set (chrom, start, end) with the cpgs, and `cov`, `meth`, `unmeth` and `avg`
#' additional fields
#' 
#' @export
smat.cpg_avg_marginals <- function(smat, cols=NULL){
    ids <- cols2ids(smat, cols)

    cov_mars <- rowSums(smat$cov[, ids])
    meth_mars <- rowSums(smat$meth[, ids])

    res <- smat$intervs %>%
        mutate(cov = cov_mars, meth = meth_mars, unmeth = cov - meth, avg = meth / cov)
    return(res)
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

    ntot <- combn(smat.colnames(smat)[ids], 2) %>% t() %>% tbl_df %>% purrr::set_names(c('cell1', 'cell2')) %>% inner_join(ntot, by=c('cell1', 'cell2'))
    
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

    rownames(smat$cov) <- NULL

    ntot <- tcrossprod(smat$cov[cgs, ids]) %>% broom::tidy()
    
    # add explicit zeroes and remove double counting
    res <-  combn(1:length(cgs), 2) %>%
        t() %>% 
        tbl_df() %>%
        purrr::set_names(c('row', 'column')) %>% 
        left_join(ntot, by=c('row', 'column')) %>% 
        mutate(value = if_else(is.na(value), 0, value))
    res <- res %>% 
        count(value) %>% 
        rename(ncells=value, npairs=n)

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
        new_mat_intervs <- smat$intervs %>% filter(id %in% ids)
    } else if (!is.null(cpg_intervs)) {
        new_mat_intervs <- smat$intervs %>% 
            inner_join(cpg_intervs, by=c('chrom', 'start', 'end')) %>% 
            arrange(id)
        ids <- new_mat_intervs$id
    } else {
        new_mat_intervs <- smat$intervs
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
        plyr::ddply(plyr::.(group), summarise_per_group, .parallel=parallel) %>%
        tbl_df %>%
        filter(ncpgs > 0) %>%
        purrr::set_names(c(group_name, 'cell', 'ncpgs', 'meth', 'unmeth', 'avg'))

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

has_stats <- function(smat){
    has_name(smat, 'stats')    
}

has_cell_metadata <- function(smat){
    has_name(smat, 'cell_metadata')    
}

has_cpg_metadata <- function(smat){
    has_name(smat, 'cpg_metadata')    
}


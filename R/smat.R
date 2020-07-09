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
    tmeth <- broom::tidy(smat$meth) %>% as_tibble() %>% rename(id=row, cell=column, meth=value)
    tunmeth <- broom::tidy(smat$unmeth) %>% as_tibble() %>% rename(id=row, cell=column, unmeth=value)
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
    df <- bind_rows(df_l, df_r) %>% inner_join(intervs, by=c('chrom', 'start', 'end')) %>% arrange(id)
    
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
#' @param smat_list list of smat objects
#' @inheritParams smat.join
#' 
#' @return joined smat object
#' 
#' @export
smat.multi_join <- function(smat_list, type='full'){
    join_func <- switch(type, 
        inner = gintervals.intersect, 
        full = gintervals.union, 
        left = function(l, r) gintervals.intersect(gintervals.union(l,r), l), 
        right = function(l, r) gintervals.intersect(gintervals.union(l,r), r))
    message(sprintf('merging intervals, join type: %s', type))    
    intervs <- purrr::reduce(smat_list, function(.x, .y) join_func(.x, .y$intervs), .init=smat_list[[1]]$intervs)
    intervs <- intervs %>% mutate(id = 1:n()) %>% as_tibble()

    message('merging matrices')
    df <- purrr::map_dfr(smat_list, ~ smat.to_df(.x) %>% select(-id))
    df <- df %>% inner_join(intervs, by=c('chrom', 'start', 'end')) %>% arrange(id)

    message('creating smat from df')
    smat <- smat.from_df(df, intervs=intervs)

    smat$stats <- purrr::map_dfr(smat_list, 'stats')
    smat$cell_metadata <- purrr::map_dfr(smat_list, 'cell_metadata')
    return(smat)
}

#' join multiple smat objects from path
#' 
#' @param paths path prefixes of the location of the smats files
#' @inheritParams smat.join
#' 
#' @return joined smat object
#' 
#' @export
smat.multi_join_from_files <- function(paths, type='full'){
    smat_list <- map(paths, smat.load)
    smat.multi_join(smat_list, type=type)
}

#' Save smat object to disk
#'
#' @param smat smat object
#'
#' @param prefix path prefix of the location of the smat files
#' @param create_dirs recursively create directories in \code{prefix}
#'
#' @return inivisibly returns the smat object
#'
#' @export
smat.save <- function(smat, prefix, create_dirs=TRUE){
    if (create_dirs){
        system(glue('mkdir -p {dirname(prefix)}'))
    }
    res <- plyr::alply(c('meth', 'unmeth', 'cov'), 1, function(x) {
        message(sprintf('saving %s', x))
        Matrix::writeMM(smat[[x]], paste0(prefix, '.', x))
    }, .parallel=FALSE)

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

    if (has_name(smat, 'cpg_metadata')){
        message('saving cpg metadata')
        fwrite(smat$cpg_metadata, paste0(prefix, '_cpg_metadata.tsv'), sep='\t')
    }

    if (has_name(smat, 'hemi_meth')){
        message('saving hemi methylation')
        fwrite(smat$hemi_meth, paste0(prefix, '_hemi_meth.tsv'), sep='\t')
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

    cpg_metadata_file <- paste0(prefix, '_cpg_metadata.tsv')
    if (file.exists(cpg_metadata_file)){
        conf$sparse_matrix$cpg_metadata <- cpg_metadata_file
    }

    hemi_meth_file <- paste0(prefix, '_hemi_meth.tsv')
    if (file.exists(hemi_meth_file)){
        conf$sparse_matrix$hemi_meth <- hemi_meth_file
    }

    tidy_cpgs_dir <- paste0(prefix, '_tcpgs')
    if (dir.exists(tidy_cpgs_dir)){
        conf$sparse_matrix$tidy_cpgs_dir <- tidy_cpgs_dir   
    }

    conf$sparse_matrix$path <- prefix

    return(smat.from_conf(conf))
}

#' Fast load sparse matrix from MatrixMarket format using data.table::fread
#'
#' @param prefix path prefix of the location of the smat files
#' 
#' @return smat object
#' 
#' @inheritDotParams Matrix::sparseMatrix
#'
#' @export
fast_readMM <- function(fn, ...){
    mm <- fread(fn, skip=1, header=FALSE)

    if (ncol(mm) == 2){
        sm <- Matrix::sparseMatrix(i = mm[, 1], j=mm[, 2])    
    } else {
        sm <- Matrix::sparseMatrix(i = mm[, 1], j=mm[, 2], x = mm[, 3])
    }
    return(sm)    
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
    mat_colnames <- fread(conf$sparse_matrix$colnames) %>% as_tibble() %>% pull(1)
    
    smat <- plyr::alply(c('cov', 'meth', 'unmeth'), 1, function(x) {
        loginfo('loading %s', x)
        m <- fast_readMM(conf$sparse_matrix[[x]]) * 1
        # m <- Matrix::readMM(conf$sparse_matrix[[x]]) * 1
        colnames(m) <- mat_colnames
        return(m)
          }, .parallel=FALSE)

    names(smat) <- c('cov', 'meth', 'unmeth')

    obs_cpgs <- fread(conf$sparse_matrix$intervs);
    if (!has_name(obs_cpgs, 'id')){
        obs_cpgs <- obs_cpgs %>% mutate(id = 1:n()) %>% as_tibble()
    } 
    smat$intervs <- as_tibble(obs_cpgs)

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

    if (has_name(conf$sparse_matrix, 'cpg_metadata')){
        smat$cpg_metadata <- fread(conf$sparse_matrix$cpg_metadata) %>% tbl_df()
    }

    if (has_name(conf$sparse_matrix, 'hemi_meth')){
        smat$hemi_meth <- fread(conf$sparse_matrix$hemi_meth) %>% tbl_df()
    }

    if (has_name(conf$sparse_matrix, 'tidy_cpgs_dir')){
        smat$tidy_cpgs_dir <- conf$sparse_matrix$tidy_cpgs_dir
    }

    if (has_name(conf$sparse_matrix, 'path')){
        smat$path <- conf$sparse_matrix$path
    }

    class(smat) <- 'smat'
    return(smat)
}


#' Get colnames of smat object
#' 
#' @param smat smat object
#' 
#' @export
smat.colnames <- function(smat){
    return(colnames(smat[['cov']]))
}

#' Get tidy cpgs of smat object
#' 
#' @param smat smat object
#' @param unique if TRUE - get unique tidy_cpgs 
#' 
#' @return smat with additional tidy_cpgs field, that contain tidy_cpgs with additional 'cell_id' column with the cell id
#' 
#' @export
smat.get_tidy_cpgs <- function(smat, unique=TRUE){
    if (!has_name(smat, 'tidy_cpgs_dir')){
        stop('tidy cpgs dir do not exist.')
    }        
        
    d <- tibble(cell_id = smat.colnames(smat)) %>% 
        mutate(dir = paste0(smat$tidy_cpgs_dir, '/', cell_id), exists = file.exists(dir))

    if (!all(d$exists)){
        warning(glue('some cells do not have tidy_cpgs: {paste(d$cell[!d$exists], collapse=", ")}'))
    }

    if (unique){
        tcpgs_dir <- 'tidy_cpgs_uniq'
    } else {
        tcpgs_dir <- 'tidy_cpgs'
    }

    d <- d %>% filter(exists) %>% select(-exists)

    tcpgs <- plyr::adply(d, 1, function(x) .gpatterns.get_tidy_cpgs_from_dir(glue('{x$dir}/{tcpgs_dir}'), uniq=unique), .parallel=TRUE) %>% select(-dir) %>% as_tibble()

    if (!unique){
        smat$tidy_cpgs_all <- tcpgs
    } else {
        smat$tidy_cpgs <- tcpgs    
    }    
    return(smat)
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
#' @return tibble with cells ('cell_id' field) and their marignal coverage ('cov' field)
#' @export
smat.cell_marginals <- function(smat, cols=NULL){
    ids <- cols2ids(smat, cols)
	mars <- colSums(smat$cov[, ids])

    return(tibble(cell_id=smat.colnames(smat)[ids], cov=mars))
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
            inner_join(select(cpg_intervs, chrom, start, end), by=c('chrom', 'start', 'end')) %>% 
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

    rownames(new_smat$cov) <- new_smat$intervs$id
    rownames(new_smat$meth) <- new_smat$intervs$id
    rownames(new_smat$unmeth) <- new_smat$intervs$id

    if (has_stats(smat)){
        new_smat$stats <- new_smat$stats %>% filter(cell_id %in% smat.colnames(new_smat))
    }

    if (has_cell_metadata(smat)){
        new_smat$cell_metadata <- new_smat$cell_metadata %>% filter(cell_id %in% smat.colnames(new_smat))
    }

    if (has_cpg_metadata(smat)){
        new_smat$cpg_metadata <- new_smat$cpg_metadata %>% inner_join(new_smat$intervs, by=c('chrom', 'start', 'end')) %>% select(-id)
    }

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
    # cell_filter <- between(colSums(smat$cov), min_cpgs, max_cpgs)
    # cg_filter <- between(rowSums(smat$cov), min_cells, max_cells)

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
    rownames(new_smat$cov) <- new_smat$intervs$id
    rownames(new_smat$meth) <- new_smat$intervs$id
    rownames(new_smat$unmeth) <- new_smat$intervs$id
    # rownames(new_smat$cov) <- NULL
    # rownames(new_smat$meth) <- NULL
    # rownames(new_smat$unmeth) <- NULL


    if (has_stats(smat)){
        new_smat$stats <- new_smat$stats %>% filter(cell_id %in% smat.colnames(new_smat))
    }

    if (has_cell_metadata(smat)){
        new_smat$cell_metadata <- new_smat$cell_metadata %>% filter(cell_id %in% smat.colnames(new_smat))
    }
    
    if (has_cpg_metadata(smat)){        
        new_smat$cpg_metadata <- new_smat$cpg_metadata[cg_filter, ]
    }

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
#' 'chrom', 'start', 'end', 'cell_id', 'ncpgs', 'cov', 'meth', 'unmeth', 'avg'
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
        group_by(chrom, start, end, cell_id) %>%
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
#' 'group_name', 'cell_id', 'ncpgs', 'meth', 'unmeth', 'avg'
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

        tibble(cell_id=names(cov), ncpgs = cov, meth=m, unmeth=um, avg=avg)
    }

    res <- smat$intervs %>%
        mutate(group = groups) %>%
        plyr::ddply(plyr::.(group), summarise_per_group, .parallel=parallel) %>%
        tbl_df %>%
        filter(ncpgs > 0) %>%
        purrr::set_names(c(group_name, 'cell_id', 'ncpgs', 'meth', 'unmeth', 'avg'))

    return(res)
}

#' Summarises CpGs by groups pf cells
#'
#' @param smat smat object
#'
#' @return tidy data frame with summariy statistics for each group and CpG
#'
#' @export
smat.summarise_cpgs <- function(smat, groups_df=NULL, groups=NULL, min_cells=1, max_cells=Inf){
    summarise_per_group <- function(x){
        cells <- match(x$cell_id, smat.colnames(smat))
        cells <- cells[!is.na(cells)]
        cgs <- which(between(rowSums(smat[['cov']][, cells]), min_cells, max_cells))

        tibble(id = cgs,
               m = rowSums(smat[['meth']][cgs, cells]),
               n = rowSums(smat[['cov']][cgs, cells]))
    }

    rm_group <- FALSE
    if (is.null(groups_df) || is.null(groups)){
        groups_df <- tibble(cell_id = smat.colnames(smat)) %>% mutate(group = 'group')
        groups <- 'group'
        rm_group <- TRUE
    }
    res <- groups_df %>% select(one_of('cell_id', groups)) %>% group_by_(groups) %>% do(summarise_per_group(.)) %>% ungroup %>% mutate(unmeth = n - m, avg = m / n)
    res <- res %>% left_join(smat$intervs, by='id') %>% select(chrom, start, end, id, one_of(groups), cov=n, meth = m, unmeth, avg)

    if (rm_group){
        res <- res %>% select(-group)
    }

    return(res)
}


# smat.summarise_cells <- function(smat){
#     browser()
# }

# .smat.summarise <- function(smat){
#     summarise_per_group <- function(cpgs, cells){
#         if (nrow(x) == 1){
#             m <- smat$meth[cpgs, cells]
#             um <- smat$unmeth[cpgs, cells]
#             cov <- smat$cov[cpgs, cells]
#             avg <- m / cov
#         } else {
#             m <- colSums(smat$meth[cpgs, cells])
#             um <- colSums(smat$unmeth[cpgs, cells])
#             cov <- colSums(smat$cov[cpgs, cells])
#             avg <- m / cov
#     }
# }

#' @export
smat.calc_cpg_metadata <- function(smat, groot=NULL, promoters_upstream=500, promoters_downstream=50, tracks=NULL, names=NULL){
    if (!is.null(groot)){
        gsetroot(groot)
    }    
    gopt <- getOption('gmultitasking')
    gdopt <- getOption('gmax.data.size')
    on.exit(options(gmultitasking = gopt))
    on.exit(options(gmax.data.size = gdopt))
    options(gmax.data.size=1e9)
    options(gmultitasking = FALSE)

    message('extracting promoters...')
    promoters <- gpatterns.get_promoters(upstream=promoters_upstream, downstream=promoters_downstream)
    gvtrack.create('prom_dist', promoters, func='distance')
    message('extracting cpg content...')
    if (!is.null(tracks)){
        if(is.null(names)) {
		names = tracks
	}
        tracks <- c(gpatterns:::.gpatterns.cg_cont_500_track, 'prom_dist', tracks)        
        names <- c('cg_cont', 'prom_dist', names)
    } else {
        tracks <- c(gpatterns:::.gpatterns.cg_cont_500_track, 'prom_dist')
        names <- c('cg_cont', 'prom_dist')
    }
    cpg_metadata <- gextract(tracks, intrevals=smat$intervs, iterator=smat$intervs, colnames=names) %>% arrange(intervalID)
    smat$cpg_metadata <- cpg_metadata %>% select(-intervalID) %>% mutate(promoter = prom_dist == 0) %>% as_tibble()
    if (has_name(smat, 'cov')){
        smat$cpg_metadata$cov <- rowSums(smat$cov)    
    }
    if (has_name(smat, 'meth')){
        smat$cpg_metadata$meth <- rowSums(smat$meth)
    }
    if (has_name(smat, 'unmeth')){
        smat$cpg_metadata$unmeth <- rowSums(smat$unmeth)
    }
    if (has_name(smat, 'cov') && has_name(smat, 'meth')){
        smat$cpg_metadata <- smat$cpg_metadata %>% mutate(avg = meth / cov)
    }

    return(smat)    
}


# smat.calc_cpg_metadata <- memoise::memoise(.calc_cpg_metadata)

smat.calc_cell_metadata <- function(smat){
    smat$cell_metadata <- tibble(cell_id = colnames(smat))
    return(smat)
}

#' @export
smat.add_intervals <- function(smat, intervals){
    # intervals <- gintervals.union(smat$intervs, intervals) %>% arrange(chrom, start, end) %>% as_tibble()
    smat_df <- smat.to_df(smat) %>% inner_join(intervals, by=c('chrom', 'start', 'end'))
    smat_new <- smat.from_df(smat_df %>% bind_rows(intervals %>% mutate(cell = 'dummy', meth = NA, unmeth=NA, cov=NA)))
    smat_new <- smat.filter(smat_new, cols=colnames(smat_new)[colnames(smat_new) != 'dummy'])
    return(smat_new)
}

#' @export
get_hemi_meth <- function(tcpgs){
    calc_hemi_meth <- function(x){
        id <- x$cell_id[1]
        message(glue('doing {id}'))        
        x <- x %>% group_by(cell_id, cg_pos) %>% filter(n_distinct(strand) >= 2)
        if (nrow(x) == 0){
            return(tibble(cell_id = id, chrom = character(), cg_pos=numeric(), `-`=numeric(), `+`=numeric()))
        }
        x <- x %>% group_by(cell_id, chrom, cg_pos, strand) %>% summarise(meth = sum(meth), cov = n()) %>% ungroup() %>% filter(cov == 1) %>% select(-cov)
        if (nrow(x) == 0){
            return(tibble(cell_id = id, chrom = character(), cg_pos=numeric(), `-`=numeric(), `+`=numeric()))   
        }
        x %>% spread(strand, meth)
    }
    res <- tcpgs %>% plyr::dlply(plyr::.(cell_id), calc_hemi_meth, .parallel=TRUE)    
    res <- res %>% map_dfr(~ .x) %>% rename(start = 'cg_pos', minus='-', plus='+') %>% mutate(end = start + 1) %>% select(cell_id, chrom, start, end, minus, plus)
    return(res)
}

#' @export
smat.calc_hemi_meth <- function(smat){
    if (!has_name(smat, 'tidy_cpgs')){
        browser()
        smat <- smat.get_tidy_cpgs(smat)
    }
    smat$hemi_meth <- get_hemi_meth(smat$tidy_cpgs)
    return(smat)
}


# Boolean functions
has_stats <- function(smat){
    has_name(smat, 'stats')  && !is.null(smat$stats)    
}

has_cell_metadata <- function(smat){
    has_name(smat, 'cell_metadata')  && !is.null(smat$cell_metadata)  
}


has_cpg_metadata <- function(smat){
    has_name(smat, 'cpg_metadata')  && !is.null(smat$cpg_metadata)   
}


#' Print smat object
#' 
#' @param smat smat object
#' 
#' @export
smat.info <- function(smat){
    ncells <- comify(ncol(smat[['cov']]))
    ncpgs <- comify(nrow(smat[['cov']]))    
    message(glue('smat object\n{ncpgs} CpGs X {ncells} cells'))    
    message(glue('--- name: {smat$name}'))    
    message(glue('--- description: "{smat$description}"'))    

    if (has_cell_metadata(smat)) {
        cell_metadata_fields <- paste(colnames(smat$cell_metadata), collapse = ', ')        
        message(glue('--- cell_metadata with the following fields: {cell_metadata_fields}'))       
        cell_groups <- group_vars(smat$cell_metadata)
        if (length(cell_groups) > 0){
            message(glue('------- cell groups: {paste(cell_groups, collapse = ", ")}')) 
        }
    }

    if (has_cpg_metadata(smat)) {
        cpg_metadata_fields <- paste(colnames(smat$cpg_metadata), collapse = ', ')        
        message(glue('--- cpg_metadata with the following fields: {cpg_metadata_fields}'))           
        cpg_groups <- group_vars(smat$cpg_metadata)
        if (length(cpg_groups) > 0){
            message(glue('------- cpg groups: {paste(cpg_groups, collapse = ", ")}')) 
        }
    }
    attrs <- paste(names(smat)[!(names(smat) %in% c('cell_metadata', 'cpg_metadata'))], collapse=', ')    
    message(glue('--- other attributes: {attrs}'))   
}

# Generics
print.smat <- function(x) smat.info(x)
colnames <-  function(x, ...) UseMethod("colnames")
colnames.default <- base::colnames
colnames.smat <- function(x) smat.colnames(x)
select.smat <- function(x, ...) smat.select(x, ...)


check_cgdb <- function(object){
    errors <- character()
    if (!all(has_name(object@cpgs, c('chrom', 'start', 'end', 'id')))){
        errors <- c(errors, 'cpgs data frame should have "chrom", "start", "end" and "id" fields')
    }

    if (!all(has_name(object@cells, c('cell_id', 'cell_num', 'plate')))){
        errors <- c(errors, 'cells data frame should have "cell_id", "cell_num" and "plate" fields')
    }

    if (!dir.exists(object@db_root)){
        errors <- c(errors, glue('root dir {object@db_root} doesn\'t exists'))
    }

    if (!dir.exists(file.path(object@db_root, 'data'))){
        errors <- c(errors, glue('data dir {object@db_root} doesn\'t exists'))
    }

    if (length(errors) == 0) TRUE else errors
}

cgdb_validate_cells <- function(db){
    db_cells <- list.files(file.path(db@db_root, 'data') , recursive=TRUE, pattern='*.idx.bin$') %>% gsub('\\.idx\\.bin', '', .) %>% gsub('/', '.', .)

    if (any(!(db@cells$cell_id %in% db_cells))){
        cells <- db@cells$cell_id[!(db@cells$cell_id %in% db_cells)]
        stop(glue('the following {length(cells)} cells do not exist in the db: {paste(cells, collapse=", ")}'))    
    }

    if (any(!(db_cells %in% db@cells$cell_id))){
        cells <- db_cells[!(db_cells %in% db@cells$cell_id)]
        warning(glue('the following {length(cells)} cells exist in the db but do not exist in db@cells: {paste(cells, collapse=", ")}'))     
    }
    invisible(db_cells)
}


#' cgdb
#'
#' @slot db_root 
#' @slot cpgs 
#' @slot cells 
#' @slot plates
#'
#' @exportClass cgdb
#' 
cgdb <- setClass(
    "cgdb",
    slots = c(
        db_root = "character",
        cpgs = "tbl_df",
        cells = "tbl_df",
        CPG_NUM = "numeric",
        .xptr = "externalptr"),
        validity = check_cgdb
)

setMethod(
  "initialize",
  signature = "cgdb",
  definition =
    function(.Object, db_root, cpgs, cells, CPG_NUM) {        
        .Object@db_root <- path.expand(db_root)
        .Object@cpgs <- cpgs
        .Object@cells <- cells
        .Object@CPG_NUM <- CPG_NUM
        .Object@.xptr <- .Call("CGDB__new", .Object@db_root, CPG_NUM)        

        return(.Object)
    }
)

#' Load cgdb database from disk
#' 
#' @param db_root root directory of the cgdb database
#' 
#' @export
cgdb_load <- function(db_root){
    cpgs_file <- glue('{db_root}/cpgs.csv')
    cells_file <- glue('{db_root}/cells.csv')
    if (!file.exists(cpgs_file)){
        stop(glue('CpGs file (cpgs.csv) doesn\'t exist. To create a new database, please run cgdb_init("{db_root}")'))
    }
    if (!file.exists(cells_file)){
        stop(glue('cells file (cells.csv) doesn\'t exist. To create a new database, please run cgdb_init("{db_root}")'))
    }
    cpgs <- fread(cpgs_file) %>% as_tibble()
    cells <- fread(cells_file) %>% as_tibble()
    db <- new('cgdb', db_root = db_root, cpgs = cpgs, cells=cells, CPG_NUM=max(cpgs$id))    
    
    return(db)    
}

#' Write cgdb to disk
#' 
#' @param db cgdb object
#' @param path path of the cgdb view
#' @param force force overwrite
#' 
#' @export
cgdb_save <- function(db, path=NULL, force=FALSE){    
    if (is.null(path)){
        path <- db@db_root
    } else {
        dir.create(path, showWarnings = FALSE, recursive = TRUE)
    }
    
    l <- flock::lock(glue('{path}/.cells_lock'))

    cells_file <- glue('{path}/cells.csv')

    if (file.exists(cells_file)){
        if (!force){
            response <- readline(glue('Are you sure you want to overwrite "{cells_file}" (Y/N)?'))    
        } else {
            response <- 'Y'
        }    

        if (response != 'Y'){
            return(NULL)           
        } 
    }

    fwrite(db@cells, cells_file, sep=',')

    cpgs_file <- glue('{path}/cpgs.csv')

    if (file.exists(cpgs_file)){
        if (!force){
            response <- readline(glue('Are you sure you want to overwrite "{cpgs_file}" (Y/N)?'))    
        } else {
            response <- 'Y'
        }    

        if (response != 'Y'){
            return(NULL)            
        } 
    }   

    fwrite(db@cpgs, cpgs_file, sep=',')
       
    if (file.exists(file.path(path, 'data'))){
        file.remove(file.path(path, 'data'))
    }

    link <- Sys.readlink(file.path(db@db_root, 'data'))
    if (link == ''){
        link <- file.path(db@db_root, 'data')
    }

    if (link != file.path(path, 'data')){
        file.symlink(file.path(db@db_root, 'data'), file.path(path, 'data'))    
    }
    
    flock::unlock(l)   
    return(NULL) 
}

#' Initialize cgdb database
#' 
#' @param db_root root directory of the cgdb database
#' @param intervals cpg intervals (id NULL 'intervs.global.seq_CG' from the current misha database would be used)
#' @param overwrite overwrite existing db
#' 
#' @export
cgdb_init <- function(db_root, intervals=NULL, overwrite=FALSE){    
    dir.create(db_root, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(db_root, 'data'), recursive = TRUE, showWarnings = FALSE)
    cells <- tibble(cell_id = character(), cell_num = numeric(), plate = character())
    if (file.exists(glue('{db_root}/cells.csv')) && !overwrite){
        stop("db exist. Run with overwrite = TRUE to overwrite")
    }
    l <- flock::lock(glue('{db_root}/.cells_lock'))
    fwrite(cells, glue('{db_root}/cells.csv'), sep=',')

    if (is.null(intervals)){
        opt <- options(gmax.data.size=1e9)
        on.exit(options(opt))
        intervals <- gintervals.load('intervs.global.seq_CG')
    }
    
    cpgs <- intervals %>% arrange(chrom, start, end) %>% mutate(id = 1:n()) %>% as_tibble()

    if (gtrack.exists("seq.CG_500_mean")){
        gvtrack.create("cg500", "seq.CG_500_mean", "avg")    
        cpgs <- gextract.left_join('cg500', intervals=cpgs, iterator=cpgs) %>% select(-(chrom1:end1)) %>% as_tibble()
    }    

    fwrite(cpgs, glue('{db_root}/cpgs.csv'), sep=',')
    flock::unlock(l)
    db <- new('cgdb', db_root=db_root, cpgs=cpgs, cells=cells, CPG_NUM=nrow(cpgs))
    return(db)
}

cgdb_update_cells <- function(db, cells, append=FALSE){   
    if (append){
        if (!is.character(db@cells$cell_num)){
            db@cells$cell_num <- as.character(db@cells$cell_num)
        }
        if (!is.character(cells$cell_num)){
            cells$cell_num <- as.character(cells$cell_num)
        }
        cells <- bind_rows(db@cells %>% filter(!(cell_id %in% cells$cell_id)), cells)
    }
    l <- flock::lock(glue('{db@db_root}/.cells_lock'))
    fwrite(cells, glue('{db@db_root}/cells.csv'), sep=',')    
    flock::unlock(l)
    db@cells <- cells
    return(db)
}

#' Remove plate from cgdb
#' 
#' @param plate_name name of the plate
#' @param force force remove (no user prompt)
#' 
#' @export
cgdb_remove_plate <- function(db, plate_name, force=FALSE){
    response <- readline(glue('Are you sure you want to remove {plate_name} and all of it\'s cells (Y/N)?'))
    if (response == 'Y' || force){
        cells <- db@cells %>% filter(plate == plate_name) %>% pull(cell_id)
        walk(cells, ~ cgdb_remove_cell(db, .x, force=TRUE))  
        system(glue('rmdir {file.path(db@db_root, "data", plate_name)}'))
        db <- cgdb_update_cells(db, db@cells %>% filter(!(cell_id %in% cells)), append=FALSE)
    }    
    
    return(db)    
}

cgdb_remove_cell <- function(db, cell_id, force=FALSE){
    if (!force){
        response <- readline(glue('Are you sure you want to remove {cell_id} (Y/N)?'))    
    } else {
        response <- 'Y'
    }
    
    if (response == 'Y'){
        x <- stringr::str_split(cell_id, '\\.')[[1]]
        file_pref <- file.path(db@db_root, 'data', x[1], x[2])
        file.remove(glue('{file_pref}.idx.bin'))
        file.remove(glue('{file_pref}.cov.bin'))
        file.remove(glue('{file_pref}.meth.bin'))                
    }
}


cgdb_add_cell <- function(db, df, cell, plate, overwrite=TRUE){    
    dirname <- file.path(db@db_root, 'data', plate)
    fname <- file.path(dirname, cell)

    if (!file.exists(paste0(fname, '.idx.bin')) || overwrite){
        dir.create(dirname, showWarnings=FALSE, recursive=TRUE)        
        df <- df %>% left_join(db@cpgs) %>% arrange(id)
        cov_vec <- df$cov
        met_vec <- df$meth
        idxs <- df$id

        writeBin(idxs, paste0(fname, '.idx.bin'), size=4)                
        writeBin(met_vec, paste0(fname, '.meth.bin'), size=4)
        writeBin(cov_vec, paste0(fname, '.cov.bin'), size=4)
        message(glue('created {cell}'))    
    }
    return(db)
}

cgdb_add_plate_from_df <- function(db, df, cells, plate_name=NULL, overwrite=TRUE, update_cells=TRUE, verbose=TRUE){
    df <- df %>% inner_join(db@cpgs %>% select(chrom, start, end, id) ,by=c('chrom', 'start', 'end')) %>% filter(!is.na(id))

    plyr::alply(cells, 1, function(x) {            
            cell <- x$cell_id
            plate <- x$plate
            cell_num <- x$cell_num

            dirname <- file.path(db@db_root, 'data', x$plate)
            fname <- file.path(dirname, x$cell_num)

            if (!file.exists(paste0(fname, '.idx.bin')) || overwrite){
                dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

                d <- df %>% filter(cell == x$cell_id)                
                cov_vec <- d$cov
                met_vec <- d$meth
                idxs <- d$id

                writeBin(idxs, paste0(fname, '.idx.bin'), size=4)                
                writeBin(met_vec, paste0(fname, '.meth.bin'), size=4)
                writeBin(cov_vec, paste0(fname, '.cov.bin'), size=4)
                if (verbose){
                    message(glue('created {cell}'))     
                }
                
            } 
            
        }, .parallel=TRUE)
    
    if (update_cells){
        db@cells <- fread(glue('{db@db_root}/cells.csv')) %>% as_tibble()
        db <- cgdb_update_cells(db, cells, append=TRUE)    
    } else {
        warning('cells.csv not updated')
    }    
    return(db)

}

#' Add sc5mc data of a plate from smat object
#' 
#' @param db cgdb object
#' @param smat smat object
#' @param plate_name name of the plate
#' @param overwrite overwrite
#' @param update_cells update cells.csv
#' @param verbose verbose messages
#' 
#' @export
cgdb_add_plate <- function(db, smat, plate_name=NULL, overwrite=TRUE, update_cells=TRUE, verbose=TRUE){
    if (is.character(smat)){
        smat <- smat.load(smat)
    }
    if (is.null(plate_name)){
        if (is.null(smat$name)){
            stop('Please provide a plate name') 
        } else {
            plate_name <- smat$name
        }       
    }

    plate_fn <- glue('{db@db_root}/{plate_name}')
    
    cells <- smat$cell_metadata %>% 
        left_join(tibble(cell_id = colnames(smat)), by = "cell_id") %>% 
        select(-one_of('plate')) %>% 
        separate(cell_id, c('plate', 'cell_num'), sep='\\.', remove=FALSE) %>%          
        # mutate(cell_num = as.numeric(cell_num)) %>% 
        arrange(cell_num)

    if (has_stats(smat)){        
        cells <- cells %>% left_join(smat$stats)
    }

    message('Converting smat to data frame')
    smat_df <- smat.to_df(smat) %>% select(-id)

    db <- cgdb_add_plate_from_df(db, smat_df, cells, plate_name=plate_name, overwrite=overwrite, update_cells=update_cells, verbose=verbose)
   
    return(db)
}

#' @export
freemem <- function(db){
    freemem_cpp(db@.xptr)
}
    
extract <- function(.Object, cells=NULL, cpgs=NULL, tidy=FALSE) {
    UseMethod("extract")
}

#' Extract data from intervals and cells
#'
#' @param cells cells to extract
#' @param cpgs cpgs to extract
#' @param tidy return tidy outout
#'
#' @export
extract.cgdb <- function(.Object, cells=NULL, cpgs=NULL, tidy=FALSE){
    if (is.null(cpgs)){
        cpgs <- .Object@cpgs    
    }   

    if (is.null(cells)){
        cells <- .Object@cells$cell_id
    }
    
    res <- extract_sc_data(.Object@.xptr, cpgs$id, cells)

    res$cov <- bind_cols(cpgs %>% select(chrom, start, end), res$cov %>% as_tibble() %>% set_names(cells) )
    res$meth <- bind_cols(cpgs %>% select(chrom, start, end), res$meth %>% as_tibble() %>% set_names(cells) )
    
    if (tidy){
        # in the future - call extract sparse
        # res <- extract_sc_data_sparse(.Object@.xptr, cpgs$id, cells)
        res$cov <- res$cov %>% gather('cell_id', 'cov', -(chrom:end))
        res$meth <- res$meth %>% gather('cell_id', 'meth', -(chrom:end))
        res <- res$cov %>% mutate(meth = res$meth$meth) %>% filter(cov > 0) 
    }
    return(res)  
}

#' Extract data from intervals and cells in a sparse format
#'
#' @param cells cells to extract
#' @param cpgs cpgs to extract
#'
#' @export
extract_sparse <- function(.Object, cells=NULL, cpgs=NULL){
    if (is.null(cpgs)){
        cpgs <- .Object@cpgs    
    }   

    if (is.null(cells)){
        cells <- .Object@cells$cell_id
    }

    res <- extract_sc_data_sparse(.Object@.xptr, cpgs$id, cells)
    res <- res %>% left_join(cpgs %>% select(chrom, start, end, id)) %>% select(chrom, start, end, cell_id, cov, meth) %>% as_tibble()
    return(res)
}

#' Summarise data from intervals and cells
#'
#' @export
summarise.cgdb <- function(.Object, tidy=TRUE){
    has_cell_groups <- is_grouped_df(.Object@cells)
    has_cpg_groups <- is_grouped_df(.Object@cpgs)

    # summary of groups of cells per CpG
    if (has_cell_groups && !has_cpg_groups){        
        res <- summarise_cells_by_group(.Object, tidy=tidy) 
    }

    # summary of groups of CpGs per cell
    if (!has_cell_groups && has_cpg_groups){   
        res <- summarise_by_cpg_groups(.Object, tidy=tidy)
    }

    # summary of groups of CpGs and cells
    if (has_cell_groups && has_cpg_groups){
        res <- summarise_by_both_groups(.Object)        
    }

    # summary of CpGs per cell
    if (!has_cell_groups && !has_cpg_groups){        
        warning('Summarising CpGs. If you want to summarise cells please call summarise_cells explicitly')
        res <- summarise_cpgs(.Object, .Object@cpgs) 
    }        

    if (is.data.frame(res)){
        res <- res %>% ungroup()  
        res <- as_tibble(res)  
    }
    
    return(res)
}


#' Summary of CpGs per cell
#' 
#' @export
summarise_cpgs <- function(db, cpgs=NULL){          
    cpgs <- cpgs %||% db@cpgs
    scdata <- mean_meth(db@.xptr, cpgs$id, db@cells$cell_id) %>% as_tibble() %>% rename(cell_id = cell)        
    if (is_grouped_df(db@cells)){
        scdata <- db@cells %>% left_join(scdata, by='cell_id')  %>% summarise(cov = sum(cov, na.rm=TRUE), meth = sum(meth, na.rm=TRUE))
    }

    return(scdata)
}

summarise_cells_by_group <- function(db, tidy=TRUE, add_metadata=TRUE){
    if (tidy){
        res <- db@cells %>% do(summarise_cells(db, .$cell_id, add_metadata=add_metadata)) %>% ungroup() %>% select(chrom, start, end, everything())  
    } else {        
        d <- db@cells %>% ungroup() %>% plyr::dlply(cell_groups(db), function(x) summarise_cells(ungroup_cells(db), x$cell_id, add_metadata=FALSE), .parallel = TRUE)
        res <- list()
        res$cov <- do.call('cbind', map(d, ~ .x[[2]])) %>% as.matrix()
        res$meth <- do.call('cbind', map(d, ~ .x[[3]])) %>% as.matrix()                
        res$intervs <- db@cpgs        
    }

    return(res)
}

#' Summary of cells per CpG
#' 
#' @export
summarise_cells <- function(db, cells=NULL, add_metadata=TRUE){
    if (is_grouped_df(db@cpgs)){
        return(as_tibble(summarise_by_both_groups(db)))
    } else {
        cells <- cells %||% db@cells$cell_id    
        scdata <- mean_meth_per_cpg(db, idxs=db@cpgs$id, cells=cells)   
        if (add_metadata)     {
            scdata <- db@cpgs %>% select(-id) %>% bind_cols(scdata)
        }        
    }
    
    return(as_tibble(scdata))
}

summarise_by_cpg_groups <- function(db, tidy=TRUE){        
    bins <- db@cpgs %>% group_indices()    
    
    res <- bin_meth_per_cell(db, db@cpgs$id, bins, db@cells$cell_id, tidy=tidy)     
    
    cpgs <- db@cpgs
    cpgs$bin <- bins
    if (tidy){
        res <- cpgs %>% distinct(bin) %>% right_join(res, by='bin') %>% select(-bin) %>% select(cell, everything()) %>% ungroup() %>% rename(cell_id = cell)   
    } else {
        res$cov <-  cpgs %>% distinct(bin) %>% right_join(res$cov, by='bin') %>% select(-bin) %>% ungroup()  %>% as_tibble()
        res$meth <-  cpgs %>% distinct(bin) %>% right_join(res$meth, by='bin') %>% select(-bin) %>% ungroup()  %>% as_tibble()        
    }
    
    return(res)  
}

summarise_by_both_groups <- function(db){
    bins <- db@cpgs %>% group_indices()    

    res <- db@cells %>% do(bin_meth(db@.xptr, db@cpgs$id, bins, .$cell_id)) %>% ungroup()

    cpgs <- db@cpgs
    cpgs$bin <- bins
    
    res <- cpgs %>% distinct(bin) %>% right_join(res, by='bin') %>% select(-bin) %>% ungroup() 
    return(res)
}

#' Summarise for intervals set
#' 
#' @param db cgdb object
#' @param intervals intervals set to summarise the genome with. Numeric values would calculate
#' the summary for equally sized bins. 
#' 
#' @return coverage and methylated calls for the intervals set. 
#' 
#' @export
summarise_intervals <- function(db, intervals){
    if (is.numeric(intervals)){
        intervals <- giterator.intervals(iterator=intervals)
    }

    res <- db %>% 
        gintervals.neighbors_cpgs(intervals) %>% 
        filter_cpgs(dist == 0) %>%
        group_by_cpgs(chrom1, start1, end1) %>% 
        summarise()
    res <- res %>% rename(chrom = chrom1, start = start1, end = end1)
    return(res)
}


#' Calculate number of cells per interval
#' 
#' @param db cgdb object
#' @param intervals intervals set to summarise the genome with. Numeric values would calculate
#' the number of cells for equally sized bins. 
#' @param min_cov minimal coverage per cell
#' 
#' @return number of cells that were covered (coverage >= \code{min_cov}) for the intervals set. 
#' 
#' @export
intervals_cell_cov <- function(db, intervals, min_cov=1){
    if (is.numeric(intervals)){
        intervals <- giterator.intervals(iterator=intervals)
    }    

    bin_intervs <- db %>% gintervals.neighbors_cpgs(intervals) %>% .@cpgs %>% group_by(chrom1, start1, end1)  
    bins <- bin_intervs %>% group_indices() 

    bin_intervs <- bin_intervs %>% ungroup() %>% mutate(bin = bins) %>% distinct(chrom1, start1, end1, bin) %>% select(chrom=chrom1, start=start1, end=end1)
   
    cov_tab <- bin_meth_per_cell_cpp(db@.xptr, db@cpgs$id, bins, db@cells$cell_id)$cov      
    colnames(cov_tab) <- db@cells$cell_id
    cov_tab <- cov_tab >= min_cov

    # in case we have only one cell in a group
    row_sums <- function(x, ...) {
        if (is.null(dim(x))) {
            return(as.numeric(x))
        }
        if (length(dim(x)) < 2){
            return(as.numeric(x))
        }      
        
        return(rowSums(x, ...))        
    }

    res <- db@cells %>% do(mutate(bin_intervs, cell_cov = row_sums(cov_tab[, .$cell_id])[-1])) %>% ungroup()
    return(res)
}

# CPP helper functions
mean_meth_per_cpg <- function(db, idxs, cells){
    bins <- 1:length(idxs)       
    res <- bin_meth(db@.xptr, idxs, bins, cells)        
    res <- res %>% mutate(id = idxs) %>% select(id, cov, meth)
    return(res)
}

bin_meth_per_cell <- function(db, idxs, bins, cells, tidy=TRUE){    
    res <- bin_meth_per_cell_cpp(db@.xptr, idxs, bins, cells)
    res$cov <- res$cov[-1, ]
    res$meth <- res$meth[-1, ]
    
    if (tidy){
        res <- purrr::map2(res[1:2], names(res[1:2]), ~ 
            as.data.frame(t(.x), row.names=cells) %>% 
            rownames_to_column() %>% 
            set_names(c('cell', res$bin)) %>% 
            gather('bin', 'var', -cell) %>% 
            set_names(c('cell', 'bin', .y)) %>%
            as_tibble())
        res <- res$cov %>% mutate(meth = res$meth$meth, bin = as.numeric(bin))    
    } else {
        res$cov <- as.data.frame(res$cov) %>% set_names(cells) %>% mutate(bin = res$bin) %>% select(bin, everything())
        res$meth <- as.data.frame(res$meth) %>% set_names(cells) %>% mutate(bin = res$bin) %>% select(bin, everything())        
    }
    
    return(res)    
}


# dplyr like functions
# group_by
#' @export
group_by_cpgs <- function(db, ...){
    db@cpgs <- db@cpgs %>% group_by(...)
    return(db)
}
#' @export
group_by_cells <- function(db, ...){
    db@cells <- db@cells %>% group_by(...)
    return(db)
}

# ungroup
#' @export
ungroup_cpgs <- function(db, ...){
    db@cpgs <- db@cpgs %>% ungroup(...)
    return(db)
}

#' @export
ungroup_cells <- function(db, ...){
    db@cells <- db@cells %>% ungroup(...)    
    return(db)
}

#' @export
ungroup.cgdb <- function(x, ...){
    x %>% ungroup_cells() %>% ungroup_cpgs()
}

# get groups
#' @export
cpg_groups <- function(db){    
    return(group_vars(db@cpgs))    
}

#' @export
cell_groups <- function(db){    
    return(group_vars(db@cells))        
}

# mutate
#' @export
mutate_cpgs <- function(db, ...){
    db@cpgs <- db@cpgs %>% mutate(...)
    return(db)
}

#' @export
g_mutate_cpgs <- function(db, track_exprs, colnames=NULL, ...){
    gdata <- gextract(track_exprs, iterator=db@cpgs, intervals=db@cpgs, colnames=colnames) %>% arrange(intervalID) %>% select(-intervalID, -(chrom:end))    
    db@cpgs <- bind_cols(db@cpgs, gdata)
    return(db)
}

#' @export
mutate_cells <- function(db, ...){
    db@cells <- db@cells %>% mutate(...)
    return(db)
}

gintervals.neighbors_cpgs <- function(db, ...){
    opt <- options(gmax.data.size = 1e9)
    on.exit(options(opt))
    db@cpgs <- db@cpgs %>% gintervals.neighbors1(...)
    return(db)
}

# filter
#' @export
filter_cpgs <- function(db, ...){
    db@cpgs <- db@cpgs %>% filter(...)    
    return(db)
}

#' @export
filter_cpgs_by_cov <- function(db, min_cov=0, max_cov=Inf){
    if (min_cov > 0 || max_cov < Inf){
        cpg_mars <- db %>% ungroup() %>% summarise_cells()
        ids <- cpg_mars %>% filter(cov >= min_cov, cov <= max_cov) %>% pull(id)        
        db <- db %>% filter_cpgs(id %in% ids)
    }
    return(db)    
}

#' Filter cgdb object by cpgs / cells coverage
#'
#' @param db smat object
#'
#' @param min_cpgs minimal number of CpGs
#' @param max_cpgs maximal number of CpGs
#' @param min_cells minimal number of cells
#' @param max_cells maximal number of cells
#'
#' @export
filter_by_cov <- function(db, min_cpgs=0, max_cpgs=Inf, min_cells=0, max_cells=Inf){
    if (min_cells > 0 || max_cells < Inf){
        cpg_mars <- db %>% ungroup() %>% summarise_cells()
        ids <- cpg_mars %>% filter(cov >= min_cells, cov <= max_cells) %>% pull(id)  
    } else {
        ids <- db@cpgs$id
    }

    if (min_cpgs > 0 || max_cpgs < Inf){
        cell_mars <- db %>% ungroup() %>% summarise_cpgs()
        cells <- cell_mars %>% filter(cov >= min_cpgs, cov <= max_cpgs) %>% pull(cell_id)
    } else {
        cells <- db@cells$cell_id
    }

    db <- db %>% filter_cells(cell_id %in% cells) %>% filter_cpgs(id %in% ids)
    return(db)
}

#' Filter cgdb object by cpgs average methylation
#'
#' @param db smat object
#'
#' @param min_avg minimal average methylation per CpG
#' @param max_avg maximal average methylation per CpG
#' @param min_cov minimal coverage per CpG
#' @param max_cov maximal coverage per CpG
#'
#' @export
filter_by_avg <- function(db, min_avg=0, max_avg=1, min_cov=1, max_cov=Inf){
    if (min_avg > 0 || max_avg < 1){
        cpg_mars <- db %>% ungroup() %>% summarise_cells()
        ids <- cpg_mars %>% 
            mutate(avg = meth / cov) %>% 
            filter(cov >= min_cov, cov <= max_cov, avg >= min_avg, avg <= max_avg) %>% 
            pull(id)  
    } else {
        ids <- db@cpgs$id
    }

    db <- db %>% filter_cpgs(id %in% ids)
    return(db)    
}

#' Filter CpGs by genomic intervals
#' 
#' @inheritDotParams gpatterns::gintervals.filter
#' @export
filter_intervals <- function(db, intervals, ...){
    opt <- options(gmax.data.size=1e9)
    on.exit(options(opt))    
    if (is.data.frame(intervals)){
        intervals <- ungroup(intervals)    
    }    
    db@cpgs <- db@cpgs %>% gpatterns::gintervals.filter(intervals)    
    return(db)
}

#' @export
filter_cells <- function(db, ...){
    db@cells <- db@cells %>% filter(...)    
    return(db)
}

#' @export
filter_cells_by_cov <- function(db, min_cov=0, max_cov=Inf){
    if (min_cov > 0 || max_cov < Inf){
        cell_mars <- db %>% ungroup() %>% summarise_cpgs()
        cells <- cell_mars %>% filter(cov >= min_cov, cov <= max_cov) %>% pull(cell_id)
        db <- db %>% filter_cells(cell_id %in% cells)        
    }
    return(db)    
}


# slice
#' @export
slice_cpgs <- function(db, ...){
    db@cpgs <- db@cpgs %>% slice(...)    
    return(db)
}


#' @export
slice_cells <- function(db, ...){
    db@cells <- db@cells %>% slice(...)    
    return(db)
}

# select
#' @export
select_cpgs <- function(db, ...){
    db@cpgs <- db@cpgs %>% select(...)    
    return(db)
}


#' @export
select_cells <- function(db, ...){
    db@cells <- db@cells %>% select(...)    
    return(db)
}

# rename
#' @export
rename_cpgs <- function(db, ...){
    db@cpgs <- db@cpgs %>% rename(...)    
    return(db)
}

#' @export
rename_cells <- function(db, ...){
    db@cells <- db@cells %>% rename(...)    
    return(db)
}


# left_join
#' @export
left_join_cpgs <- function(db, ...){
    db@cpgs <- db@cpgs %>% left_join(...)    
    return(db)
}

#' @export
left_join_cells <- function(db, ...){
    db@cells <- db@cells %>% left_join(...)    
    return(db)
}

#' @export
inner_join_cpgs <- function(db, ...){
    db@cpgs <- db@cpgs %>% inner_join(...)    
    return(db)
}

#' @export
inner_join_cells <- function(db, ...){
    db@cells <- db@cells %>% inner_join(...)    
    return(db)
}

#' @export
anti_join_cpgs <- function(db, ...){
    db@cpgs <- db@cpgs %>% anti_join(...)    
    return(db)
}

#' @export
anti_join_cells <- function(db, ...){
    db@cells <- db@cells %>% anti_join(...)    
    return(db)
}

# count
#' @export
count_cpgs <- function(db, ...){
    db@cpgs <- db@cpgs %>% count(...)    
    return(db)
}

#' @export
count_cells <- function(db, ...){
    db@cells <- db@cells %>% count(...)    
    return(db)
}

# other functions
#' @export
call_cells <- function(db, func, ...){
    db@cells <- db@cells %>% func(...)
    return(db)
}

#' @export
call_cpgs <- function(db, func, ...){
    db@cpgs <- db@cpgs %>% func(...) 
    return(db)
}

#' @export
pull_cells <- function(db){
    db@cells
}

#' @export
pull_cpgs <- function(db){
    db@cpgs
}


setMethod("show", 
          "cgdb", 
          function(object){
            cgdb_info(object)
          }
)

cgdb_info <- function(db){
    ncells <- comify(nrow(db@cells))
    ncpgs <- comify(nrow(db@cpgs))
    message(glue('cgdb object\n{ncpgs} CpGs X {ncells} cells'))    
    message(glue('--- root (@db_root): {db@db_root}'))
    if (length(cell_groups(db)) > 0){
        groups <- paste(cell_groups(db), collapse = ', ')
        message(glue('--- cell groups: {groups}'))    
    }

    if (length(cpg_groups(db)) > 0){
        groups <- paste(cpg_groups(db), collapse = ', ')
        message(glue('--- CpG groups: {groups}'))    
    }
    
    # message('\n--- Cells (@cells):')
    # print(db@cells)
    # message('\n--- CpGs (@cpgs):')
    # print(db@cpgs)
}

#' Count methylation calls (00,01,10,11) for pairs of cells
#' 
#' @param db cgdb object
#' 
#' @return data frame with the following fields: cell1, cell2, n00, n01, n10, n11
#' 
#' @export
count_pairs <- function(db, max_chunk_size=200){    
    message('filtering CpGs with less than 2 cells')
    db <- db %>% filter_by_cov(min_cells=2)

    cell_pairs <- t(combn(db@cells$cell_id, 2)) %>% as_tibble() %>% set_names(c('cell1', 'cell2'))
    cell_pairs <- arrange(cell_pairs, cell1, cell2)

    cells <- tibble(cell = db@cells$cell_id, chunk = ntile(cell, round(length(cell) / max_chunk_size)))
    
    message(glue('chunks: {length(unique(cells$chunk))}'))
    doMC::registerDoMC(min(5, getOption('gpatterns.parallel.thread_num')))
    
    pairs <- plyr::dlply(cells, plyr::.(chunk), function(x) {    
        chunk <- x$chunk[1]        
        message(glue('extracting ({chunk})'))
        df <- extract_sc_data_sparse(db@.xptr, db@cpgs$id, x$cell)    

        message(glue('arranging ({chunk})'))
        smat <- list()
        smat$meth <- tidy2smat(df, 'id', 'cell_id', 'meth')
        smat$cov <- tidy2smat(df, 'id', 'cell_id', 'cov')       
        
        colnames(smat$cov) <- x$cell
        colnames(smat$meth) <- x$cell

        print(dim(smat$meth))          

        message(glue('creating unmeth matrix ({chunk})'))        
        smat$unmeth <- smat$cov - smat$meth    

        message(glue('calculating pairs ({chunk})'))
        pairs_meth <- sc5mc.calc_pdiff(smat, min_cgs=0) %>% select(cell1, cell2, n00, n01, n10, n11)        
        rm(smat)
        gc()
        return(pairs_meth)
    }, .parallel = TRUE)

    pairs <- map_df(pairs, ~.x) 
    
    cell_pairs <- cell_pairs %>% anti_join(pairs, by=c('cell1', 'cell2')) %>% arrange(cell1, cell2)
    message(glue('remaining {nrow(cell_pairs)}'))
    
    if (nrow(cell_pairs) > 0){
        # gpatterns::gpatterns.set_parallel(getOption('gpatterns.parallel.thread_num'))
        # nbins <- getOption('gpatterns.parallel.thread_num')
        nbins <- 5

        cell_pairs <- cell_pairs %>%
            mutate(chunk = ntile(cell1, nbins)) %>% 
            group_by(cell1) %>% 
            mutate(chunk = first(chunk)) %>%
            ungroup()
        
        pairs_meth <- plyr::ddply(cell_pairs, plyr::.(chunk), function(x) count_pairs_all_cpp(db@.xptr, db@cpgs$id, x$cell1, x$cell2) %>% as_tibble() %>% set_names(c('n00', 'n01', 'n10', 'n11')), .parallel = TRUE )  %>% select(-chunk)
       
        pairs_meth <- bind_cols(cell_pairs, pairs_meth) %>% select(-chunk)
        pairs_meth <- bind_rows(pairs_meth, rename(pairs_meth, cell1 = cell2, cell2 = cell1))
        pairs_meth <- bind_rows(pairs, pairs_meth)
    } else {
        pairs_meth <- pairs
    }

    return(pairs_meth)    
}

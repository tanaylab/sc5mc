
check_cgdb <- function(object){
    errors <- character()
    if (!all(has_name(object@cpgs, c('chrom', 'start', 'end', 'id')))){
        errors <- c(errors, 'cpgs data frame should have "chrom", "start", "end" and "id" fields')
    }

    if (!all(has_name(object@cells, c('cell_id', 'cell_num', 'plate')))){
        errors <- c(errors, 'cells data frame should have "cell_id", "cell_num" and "plate" fields')
    }

    if (length(errors) == 0) TRUE else errors
}

cgdb_validate_cells <- function(cgdb){
    db_cells <- list.files(file.path(cgdb@db_root, 'data') , recursive=TRUE, pattern='*.idx.bin$') %>% gsub('\\.idx\\.bin', '', .) %>%  gsub('/', '.', .)

    if (any(!(cgdb@cells$cell_id %in% db_cells))){
        cells <- cgdb@cells$cell_id[!(cgdb@cells$cell_id %in% db_cells)]
        stop(glue('the following cells do not exist in the db: {paste(cells, collapse=", ")}'))     
    }

    if (any(!(db_cells %in% cgdb@cells$cell_id))){
        cells <- db_cells[!(db_cells %in% cgdb@cells$cell_id)]
        warning(glue('the following cells exist in the db but do not exist in cgdb@cells: {paste(cells, collapse=", ")}'))     
    }
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
        stop(glue('CpGs file (cpgs.csv) doesn\'t exist. Please run cgdb_init("{db_root}")'))
    }
    if (!file.exists(cells_file)){
        stop(glue('cells file (cells.csv) doesn\'t exist. Please run cgdb_init("{db_root}")'))
    }
    cpgs <- fread(cpgs_file) %>% as_tibble()
    cells <- fread(cells_file) %>% as_tibble()
    cgdb <- new('cgdb', db_root = db_root, cpgs = cpgs, cells=cells, CPG_NUM=nrow(cpgs))    
    
    return(cgdb)    
}

#' Write cgdb database to disk
#' 
#' @param cgdb cgdb object
#' 
#' @export
cgdb_save <- function(cgdb){
    l <- flock::lock(glue('{cgdb@db_root}/.cells_lock'))
    fwrite(cgdb@cells, glue('{cgdb@db_root}/cells.csv'), sep=',')   
    fwrite(cgdb@cpgs, glue('{cgdb@db_root}/cpgs.csv'), sep=',') 
    flock::unlock(l)
    return(cgdb)
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
    if (file.exists(glue('{db_root}/cells.csv')) || !overwrite){
        stop("db exist. Run with overwrite = TRUE to overwrite")
    }
    l <- flock::lock(glue('{cgdb@db_root}/.cells_lock'))
    fwrite(cells, glue('{db_root}/cells.csv'), sep=',')

    if (is.null(intervals)){
        opt <- options(gmax.data.size=1e9)
        on.exit(options(opt))
        intervals <- gintervals.load('intervs.global.seq_CG')
    }
    
    cpgs <- intervals %>% arrange(chrom, start, end) %>% mutate(id = 1:n()) %>% as_tibble()
    fwrite(cpgs, glue('{db_root}/cpgs.csv'), sep=',')
    flock::unlock(l)
    cgdb <- new('cgdb', db_root=db_root, cpgs=cpgs, cells=cells, CPG_NUM=nrow(cpgs))
    return(cgdb)
}

cgdb_update_cells <- function(cgdb, cells, append=FALSE){   
    if (append){
        cells <- bind_rows(cgdb@cells %>% filter(!(cell_id %in% cells$cell_id)), cells)
    }
    l <- flock::lock(glue('{cgdb@db_root}/.cells_lock'))
    fwrite(cells, glue('{cgdb@db_root}/cells.csv'), sep=',')    
    flock::unlock(l)
    cgdb@cells <- cells
    return(cgdb)
}

#' Remove plate from cgdb#' 
#' 
#' @param plate_name name of the plate
#' @param force force remove (no user prompt)
#' 
#' @export
cgdb_remove_plate <- function(cgdb, plate_name, force=FALSE){
    response <- readline(glue('Are you sure you want to remove {plate_name} and all of it\'s cells (Y/N)?'))
    if (response == 'Y' || force){
        cells <- cgdb@cells %>% filter(plate == plate_name) %>% pull(cell_id)
        walk(cells, ~ cgdb_remove_cell(cgdb, .x, force=TRUE))  
        system(glue('rmdir {file.path(cgdb@db_root, "data", plate_name)}'))
        cgdb <- cgdb_update_cells(cgdb, cgdb@cells %>% filter(!(cell_id %in% cells)), append=FALSE)
    }    
    
    return(cgdb)    
}

cgdb_remove_cell <- function(cgdb, cell_id, force=FALSE){
    if (!force){
        response <- readline(glue('Are you sure you want to remove {cell_id} (Y/N)?'))    
    } else {
        response <- 'Y'
    }
    
    if (response == 'Y'){
        x <- stringr::str_split(cell_id, '\\.')[[1]]
        file_pref <- file.path(cgdb@db_root, 'data', x[1], x[2])
        file.remove(glue('{file_pref}.idx.bin'))
        file.remove(glue('{file_pref}.cov.bin'))
        file.remove(glue('{file_pref}.meth.bin'))                
    }
}

#' Add sc5mc data of a plate from smat object
#' 
#' @param cgdb cgdb object
#' @param smat smat object
#' @param plate_name name of the plate
#' @param overwrite overwrite
#' 
#' @export
cgdb_add_plate <- function(cgdb, smat, plate_name=NULL, overwrite=TRUE){
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

    plate_fn <- glue('{cgdb@db_root}/{plate_name}')
    # if (!(all(smat$intervs$chrom == cgdb@cpgs$chrom) && all(smat$intervs$start == cgdb@cpgs$start) && all(smat$intervs$end == cgdb@cpgs$end))){                  
    #     cell_md <- smat$cell_metadata
    #     smatc<- smat.add_intervals(smat, cgdb@cpgs %>% select(chrom, start, end))   
    #     smat$cell_metadata <- cell_md
    # }    

    cells <- smat$cell_metadata %>% 
        left_join(tibble(cell_id = colnames(smat))) %>% 
        select(-one_of('plate')) %>% 
        separate(cell_id, c('plate', 'cell_num'), sep='\\.', remove=FALSE) %>%          
        # mutate(cell_num = as.numeric(cell_num)) %>% 
        arrange(cell_num)

    # stopifnot(all(smat$intervs$chrom == cgdb@cpgs$chrom) && all(smat$intervs$start == cgdb@cpgs$start) && all(smat$intervs$end == cgdb@cpgs$end))    

    smat_df <- smat.to_df(smat) %>% select(-id) %>% inner_join(cgdb@cpgs %>% select(chrom, start, end, id) ,by=c('chrom', 'start', 'end')) %>% filter(!is.na(id))    

    plyr::alply(cells, 1, function(x) {            
            cell <- x$cell_id
            plate <- x$plate
            cell_num <- x$cell_num

            dirname <- file.path(cgdb@db_root, 'data', x$plate)
            fname <- file.path(dirname, x$cell_num)

            if (!file.exists(paste0(fname, '.idx.bin')) || overwrite){
                dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

                # cov_vec <- smat$cov[, x$cell_id]
                # met_vec <- smat$meth[, x$cell_id]
                # idxs <- which(cov_vec > 0)

                d <- smat_df %>% filter(cell == x$cell_id)                
                cov_vec <- d$cov
                met_vec <- d$meth
                idxs <- d$id

                writeBin(idxs, paste0(fname, '.idx.bin'), size=4)
                # writeBin(met_vec[idxs], paste0(fname, '.meth.bin'), size=4)
                # writeBin(cov_vec[idxs], paste0(fname, '.cov.bin'), size=4)
                writeBin(met_vec, paste0(fname, '.meth.bin'), size=4)
                writeBin(cov_vec, paste0(fname, '.cov.bin'), size=4)
                message(glue('created {cell}'))    
            } 
            
        }, .parallel=FALSE)
    
    cgdb@cells <- fread(glue('{cgdb@db_root}/cells.csv')) %>% as_tibble()
    cgdb <- cgdb_update_cells(cgdb, cells, append=TRUE)
    return(cgdb)
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
        res$cov <- res$cov %>% gather('cell_id', 'cov', -(chrom:end))
        res$meth <- res$meth %>% gather('cell_id', 'meth', -(chrom:end))
        res <- res$cov %>% mutate(meth = res$meth$meth) %>% filter(cov > 0) 
    }
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
        res <- .Object@cells %>% do(summarise_cells(.Object, .$cell_id)) %>% ungroup() %>% select(chrom, start, end, everything())      
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
summarise_cpgs <- function(cgdb, cpgs=NULL){        
    cpgs <- cpgs %||% cgdb@cpgs
    scdata <- mean_meth(cgdb@.xptr, cpgs$id, cgdb@cells$cell_id) %>% as_tibble() %>% rename(cell_id = cell)    
    return(scdata)
}

#' Summary of cells per CpG
#' 
#' @export
summarise_cells <- function(cgdb, cells=NULL){
    cells <- cells %||% cgdb@cells$cell_id
    scdata <- mean_meth_per_cpg(cgdb, idxs=cgdb@cpgs$id, cells=cells)    
    scdata <- cgdb@cpgs %>% select(-id) %>% bind_cols(scdata) %>% as_tibble() 
    return(scdata)    
}

summarise_by_cpg_groups <- function(cgdb, tidy=TRUE){        
    bins <- cgdb@cpgs %>% group_indices()    
    
    res <- bin_meth_per_cell(cgdb, cgdb@cpgs$id, bins, cgdb@cells$cell_id, tidy=tidy)     
    
    cpgs <- cgdb@cpgs
    cpgs$bin <- bins
    if (tidy){
        res <- cpgs %>% distinct(bin, .by_group=TRUE) %>% right_join(res, by='bin') %>% select(-bin) %>% select(cell, everything()) %>% ungroup()    
    } else {
        res$cov <-  cpgs %>% distinct(bin, .by_group=TRUE) %>% right_join(res$cov, by='bin') %>% select(-bin) %>% ungroup()  %>% as_tibble()
        res$meth <-  cpgs %>% distinct(bin, .by_group=TRUE) %>% right_join(res$meth, by='bin') %>% select(-bin) %>% ungroup()  %>% as_tibble()        
    }
    
    return(res)  
}

summarise_by_both_groups <- function(cgdb){
    bins <- cgdb@cpgs %>% group_indices()    

    res <- cgdb@cells %>% do(bin_meth(cgdb@.xptr, cgdb@cpgs$id, bins, .$cell_id)) %>% ungroup()

    cpgs <- cgdb@cpgs
    cpgs$bin <- bins
    
    res <- cpgs %>% distinct(bin, .by_group=TRUE) %>% right_join(res, by='bin') %>% select(-bin) %>% ungroup() 
    return(res)
}

# CPP helper functions
mean_meth_per_cpg <- function(cgdb, idxs, cells){
    bins <- 1:length(idxs)        
    res <- bin_meth(cgdb@.xptr, idxs, bins, cells)
    res <- res %>% mutate(id = idxs) %>% select(id, cov, meth)
    return(res)
}

bin_meth_per_cell <- function(cgdb, idxs, bins, cells, tidy=TRUE){    
    res <- bin_meth_per_cell_cpp(cgdb@.xptr, idxs, bins, cells)
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
group_by_cpgs <- function(cgdb, ...){
    cgdb@cpgs <- cgdb@cpgs %>% group_by(...)
    return(cgdb)
}
#' @export
group_by_cells <- function(cgdb, ...){
    cgdb@cells <- cgdb@cells %>% group_by(...)
    return(cgdb)
}

# ungroup
#' @export
ungroup_cpgs <- function(cgdb, ...){
    cgdb@cpgs <- cgdb@cpgs %>% ungroup(...)
    return(cgdb)
}

#' @export
ungroup_cells <- function(cgdb, ...){
    cgdb@cells <- cgdb@cells %>% ungroup(...)    
    return(cgdb)
}

# get groups
#' @export
cpg_groups <- function(cgdb){    
    return(group_vars(cgdb@cpgs))    
}

#' @export
cell_groups <- function(cgdb){    
    return(group_vars(cgdb@cells))        
}

# mutate
#' @export
mutate_cpgs <- function(cgdb, ...){
    cgdb@cpgs <- cgdb@cpgs %>% mutate(...)
    return(cgdb)
}

#' @export
g_mutate_cpgs <- function(cgdb, track_exprs, colnames=NULL, ...){
    gdata <- gextract(track_exprs, iterator=cgdb@cpgs, intervals=cgdb@cpgs, colnames=colnames) %>% arrange(intervalID) %>% select(-intervalID, -(chrom:end))    
    cgdb@cpgs <- bind_cols(cgdb@cpgs, gdata)
    return(cgdb)
}

#' @export
mutate_cells <- function(cgdb, ...){
    cgdb@cells <- cgdb@cells %>% mutate(...)
    return(cgdb)
}

# filter
#' @export
filter_cpgs <- function(cgdb, ...){
    cgdb@cpgs <- cgdb@cpgs %>% filter(...)    
    return(cgdb)
}

#' 
#' @inheritDotParams gpatterns::gintervals.filter
#' @export
filter_intervals <- function(cgdb, intervals, ...){
    opt <- options(gmax.data.size=1e9)
    on.exit(options(opt))
    cgdb@cpgs <- cgdb@cpgs %>% gintervals.filter(intervals)    
    return(cgdb)
}

#' @export
filter_cells <- function(cgdb, ...){
    cgdb@cells <- cgdb@cells %>% filter(...)    
    return(cgdb)
}


# slice
#' @export
slice_cpgs <- function(cgdb, ...){
    cgdb@cpgs <- cgdb@cpgs %>% slice(...)    
    return(cgdb)
}


#' @export
slice_cells <- function(cgdb, ...){
    cgdb@cells <- cgdb@cells %>% slice(...)    
    return(cgdb)
}

# select
#' @export
select_cpgs <- function(cgdb, ...){
    cgdb@cpgs <- cgdb@cpgs %>% select(...)    
    return(cgdb)
}


#' @export
select_cells <- function(cgdb, ...){
    cgdb@cells <- cgdb@cells %>% select(...)    
    return(cgdb)
}



check_cgdb <- function(object){
    errors <- character()
    if (!all(has_name(object@cpgs, c('chrom', 'start', 'end', 'id')))){
        errors <- c(errors, 'cpgs data frame should have "chrom", "start", "end" and "id" fields')
    }

    if (!all(has_name(object@cells, c('cell_id', 'cell_num', 'plate')))){
        errors <- c(errors, 'cells data frame should have "cell_id", "cell_num" and "plate" fields')
    }

    # plates <- unique(object@cells$plate)
    # if (!all(file.exists(glue('{object@db_root}/{plates}.h5')))){
    #     warning('some plate files do not exist')
    # }

    if (length(errors) == 0) TRUE else errors
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
        CPG_NUM = "numeric"),
        validity = check_cgdb
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
    fwrite(cgdb@cells, glue('{cgdb@db_root}/cells.csv'), sep=',')   
    fwrite(cgdb@cpgs, glue('{cgdb@db_root}/cpgs.csv'), sep=',') 
    return(cgdb)
}

#' Initialize cgdb database
#' 
#' @param db_root root directory of the cgdb database
#' @param intervals cpg intervals (id NULL 'intervs.global.seq_CG' from the current misha database would be used)
#' 
#' @export
cgdb_init <- function(db_root, intervals=NULL){
    dir.create(db_root, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(db_root, 'data'), recursive = TRUE, showWarnings = FALSE)
    cells <- tibble(cell_id = character(), cell_num = numeric(), plate = character())
    fwrite(cells, glue('{db_root}/cells.csv'), sep=',')

    if (is.null(intervals)){
        opt <- options(gmax.data.size=1e9)
        on.exit(options(opt))
        intervals <- gintervals.load('intervs.global.seq_CG')
    }
    
    cpgs <- intervals %>% arrange(chrom, start, end) %>% mutate(id = 1:n()) %>% as_tibble()
    fwrite(cpgs, glue('{db_root}/cpgs.csv'), sep=',')

    cgdb <- new('cgdb', db_root=db_root, cpgs=cpgs, cells=cells, CPG_NUM=nrow(cpgs))
    return(cgdb)
}

cgdb_update_cells <- function(cgdb, cells, append=FALSE){   
    if (append){
        cells <- bind_rows(cgdb@cells %>% filter(!(cell_id %in% cells$cell_id)), cells)
    }
    fwrite(cells, glue('{cgdb@db_root}/cells.csv'), sep=',')    
    cgdb@cells <- cells
    return(cgdb)
}

#' Add sc5mc data of a plate from smat object
#' 
#' @param cgdb cgdb object
#' @param smat smat object
#' @param plate_name name of the plate
#' 
#' @export
cgdb_add_plate <- function(cgdb, smat, plate_name=NULL){
    if (is.null(plate_name)){
        if (is.null(smat$name)){
            stop('Please provide a plate name') 
        } else {
            plate_name <- smat$name
        }       
    }

    plate_fn <- glue('{cgdb@db_root}/{plate_name}')
    if (!(all(smat$intervs$chrom == cgdb@cpgs$chrom) && all(smat$intervs$start == cgdb@cpgs$start) && all(smat$intervs$end == cgdb@cpgs$end))){  
        print('problem')      
        browser()
        smat <- smat.add_intervals(smat, cgdb@cpgs %>% select(chrom, start, end))   
    }

    cells <- smat$cell_metadata %>% 
        left_join(tibble(cell_id = colnames(smat))) %>% 
        separate(cell_id, c('plate1', 'cell_num'), sep='\\.', remove=FALSE) %>% 
        select(-plate1) %>% 
        mutate(cell_num = as.numeric(cell_num)) %>% 
        arrange(cell_num)

    stopifnot(all(smat$intervs$chrom == cgdb@cpgs$chrom) && all(smat$intervs$start == cgdb@cpgs$start) && all(smat$intervs$end == cgdb@cpgs$end))    

    plyr::alply(cells, 1, function(x) {
            cell <- x$cell_id
            plate <- x$plate
            cell_num <- x$cell_num

            cov_vec <- smat$cov[, x$cell_id]
            met_vec <- smat$meth[, x$cell_id]

            idxs <- which(cov_vec > 0)

            dirname <- file.path(cgdb@db_root, 'data', x$plate)
            fname <- file.path(dirname, x$cell_num)

            dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

            writeBin(idxs, paste0(fname, '.idx.bin'), size=4)
            writeBin(met_vec[idxs], paste0(fname, '.meth.bin'), size=4)
            writeBin(cov_vec[idxs], paste0(fname, '.cov.bin'), size=4)
        }, .parallel=TRUE)
    

    cgdb <- cgdb_update_cells(cgdb, cells, append=TRUE)
    return(cgdb)
}
    


#' Extract data from intervals and cells
#'
#' @param cells cells to extract
#' @param cpgs cpgs to extract
#'
#' @export
extract.cgdb <- function(.Object, cells=NULL, cpgs=NULL){
    if (is.null(cpgs)){
        cpgs <- .Object@cpgs    
    }   

    if (is.null(cells)){
        cells <- .Object@cells$cell_id
    }

    return(extract_sc_data(.Object, cells, cpgs %>% select(chrom, start, end, id)))  
}

#' Summarise data from intervals and cells
#'
#' @export
summarise.cgdb <- function(.Object){
    has_cell_groups <- is_grouped_df(.Object@cells)
    has_cpg_groups <- is_grouped_df(.Object@cpgs)

    # summary of groups of cells per CpG
    if (has_cell_groups && !has_cpg_groups){        
        res <- .Object@cells %>% do(summarise_cells(.Object, .$cell_id)) %>% ungroup() %>% select(chrom, start, end, everything())      
    }

    # summary of groups of CpGs per cell
    if (!has_cell_groups && has_cpg_groups){   
        res <- summarise_by_cpg_groups(.Object)
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

    res <- res %>% ungroup()

    return(as_tibble(res))
}


#' Summary of CpGs per cell
#' 
#' @export
summarise_cpgs <- function(cgdb, cpgs=NULL){        
    cpgs <- cpgs %||% cgdb@cpgs
    scdata <- mean_meth(cpgs$id, file.path(cgdb@db_root, 'data'), cgdb@cells$cell_id, cgdb@CPG_NUM) %>% as_tibble() %>% rename(cell_id = cell)    
    return(scdata)
}

#' Summary of cells per CpG
#' 
#' @export
summarise_cells <- function(cgdb, cells=NULL){
    cells <- cells %||% cgdb@cells$cell_id
    scdata <- mean_meth_per_cpg(idxs=cgdb@cpgs$id, db_dir=file.path(cgdb@db_root, 'data'), cells=cells, CPG_NUM=cgdb@CPG_NUM)    
    scdata <- cgdb@cpgs %>% select(-id) %>% bind_cols(scdata) %>% as_tibble() 
    return(scdata)    
}

summarise_by_cpg_groups <- function(cgdb){        
    bins <- cgdb@cpgs %>% group_indices()    
    
    res <- bin_meth_per_cell(cgdb@cpgs$id, bins, file.path(cgdb@db_root, 'data'), cgdb@cells$cell_id, cgdb@CPG_NUM)     

    cpgs <- cgdb@cpgs
    cpgs$bin <- bins
    res <- cpgs %>% distinct(bin, .by_group=TRUE) %>% right_join(res, by='bin') %>% select(-bin) %>% select(cell, everything()) %>% ungroup()
    return(res)  
}

summarise_by_both_groups <- function(cgdb){
    bins <- cgdb@cpgs %>% group_indices()    

    res <- cgdb@cells %>% do(bin_meth(cgdb@cpgs$id, bins, file.path(cgdb@db_root, 'data'), .$cell_id, cgdb@CPG_NUM)) %>% ungroup()

    cpgs <- cgdb@cpgs
    cpgs$bin <- bins
    
    res <- cpgs %>% distinct(bin, .by_group=TRUE) %>% right_join(res, by='bin') %>% select(-bin) %>% ungroup() 
    return(res)
}

# CPP helper functions
mean_meth_per_cpg <- function(idxs, db_dir, cells, CPG_NUM){
    bins <- 1:length(idxs)        
    res <- bin_meth(idxs, bins, db_dir, cells, CPG_NUM)
    res <- res %>% mutate(id = idxs) %>% select(id, cov, meth)
    return(res)
}

bin_meth_per_cell <- function(idxs, bins, db_dir, cells, CPG_NUM){    
    res <- bin_meth_per_cell_cpp(idxs, bins, db_dir, cells, CPG_NUM)

    res <- purrr::map2(res[1:2], names(res[1:2]), ~ 
        as.data.frame(.x[, -1], row.names=cells) %>% 
        rownames_to_column() %>% 
        set_names(c('cell', res$bin)) %>% 
        gather('bin', 'var', -cell) %>% 
        set_names(c('cell', 'bin', .y)) %>%
        as_tibble())
    res <- res$cov %>% mutate(meth = res$meth$meth, bin = as.numeric(bin))
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

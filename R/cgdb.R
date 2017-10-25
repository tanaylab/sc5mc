
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

# extract_sc_data <- function(cgdb, cells, cpgs){
#     cells_md <- cgdb@cells %>% ungroup() %>% filter(cell_id %in% cells) %>% select(cell_id, cell_num, plate)

#     res <- cells_md %>% arrange(plate, cell_num) %>% plyr::ddply(plyr::.(plate), function(x){       
#         plate <- x$plate[1]
#         h5f <- get_h5_obj(cgdb, plate)
#         m_meth <- h5f['meth'][cpgs$id, x$cell_num]
#         m_cov <- h5f['cov'][cpgs$id, x$cell_num]
#         cpgs_f <- rowSums(m_cov) > 0

#         cov_tab <- bind_cols(cpgs[cpgs_f, ], as.data.frame(m_cov) %>% set_names(x$cell_id) %>% filter(cpgs_f)) %>% as_tibble() %>% select(-id)
#         meth_tab <- bind_cols(cpgs[cpgs_f, ], as.data.frame(m_meth) %>% set_names(x$cell_id) %>% filter(cpgs_f)) %>% as_tibble() %>% select(-id)
#         stab <- cov_tab %>% gather('cell_id', 'cov', -(chrom:end)) %>% mutate(meth = meth_tab %>% gather('cell_id', 'meth', -(chrom:end)) %>% pull(meth)) %>% filter(cov > 0)
#         return(stab)
#     } )
#     res <- res %>% select(chrom, start, end, cell_id, cov, meth) %>% as_tibble()
    
#     return(res) 
# }

#' Summarise data from intervals and cells
#'
#' @export
summarise.cgdb <- function(.Object){
    has_cell_groups <- is_grouped_df(.Object@cells)
    has_cpg_groups <- is_grouped_df(.Object@cpgs)

    if (has_cell_groups && !has_cpg_groups){        
        res <- .Object@cells %>% do(summarise_cells(.Object, .$cell_id)) %>% ungroup() %>% select(chrom, start, end, everything())      
    }
    if (!has_cell_groups && has_cpg_groups){        
        res <- .Object@cpgs %>% do(summarise_cpgs(.Object, .)) %>% ungroup()        
    }

    if (has_cell_groups && has_cpg_groups){
        f <- function(.Object, cpgs) { .Object@cells %>% do(summarise_sc_data(.Object, .$cell_id, cpgs)) }
        res <- .Object@cpgs %>% do(f(.Object, .)) %>% ungroup()     
    }

    if (!has_cell_groups && !has_cpg_groups){
        # res <- extract(.Object)
        res <- summarise_cpgs(.Object, .Object@cpgs) 
    }        

    return(res)
}

summarise_cpgs <- function(cgdb, cpgs){        
    scdata <- mean_meth(cpgs$id, file.path(cgdb@db_root, 'data'), cgdb@cells$cell_id, cgdb@CPG_NUM) %>% as_tibble() %>% rename(cell_id = cell)    
    return(scdata)
}

summarise_cells <- function(cgdb, cells){   
    browser()
    # scdata <- cgdb %>% extract(cells, cgdb@cpgs)
    # scdata <- scdata %>% group_by(chrom, start, end) %>% summarise(cov = sum(cov), meth = sum(meth)) %>% ungroup()
    return(scdata)
    
}

summarise_sc_data <- function(cgdb, cells, cpgs){
    # scdata <- cgdb %>% extract(cells, cpgs) %>% summarise(cov = sum(cov), meth = sum(meth))
    browser()
    return(scdata)
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

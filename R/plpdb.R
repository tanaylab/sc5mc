#' plpdb
#'
#' @slot db_root 
#' @slot cpgs 
#' @slot cells 
#' @slot CPG_NUM
#'
#' @exportClass cgdb
#' 
plpdb <- setClass(
    "plpdb",
    slots = c(binsize = "numeric"),
    contains = 'cgdb'
)

#' Initialize pileup database
#' 
#' @param db_root root directory of the cgdb database
#' @param binsize size of genomic bins
#' @param overwrite overwrite existing db
#' @param cgdb cgdb object to initialize from
#' 
#' @export
plpdb_init <- function(db_root, binsize, overwrite=FALSE){    
    opt <- options(gmax.data.size=1e9)    
    on.exit(options(opt))

    l <- .create_cgdb(db_root, intervals=giterator.intervals(iterator=binsize), overwrite=overwrite, add_cg_cont=FALSE)
    db <- new('plpdb', db_root=db_root, cpgs=l$cpgs, cells=l$cells, CPG_NUM=nrow(l$cpgs))
    return(db)    
}

setMethod("show", 
          "plpdb", 
          function(object){
            plpdb_info(object)
          }
)

plpdb_info <- function(db){
    ncells <- comify(nrow(db@cells))
    nintervals <- comify(nrow(db@cpgs))
    message(glue('plpdb object\n{nintervals} intervals X {ncells} cells'))    
    message(glue('--- root (@db_root): {db@db_root}'))
    if (length(cell_groups(db)) > 0){
        groups <- paste(cell_groups(db), collapse = ', ')
        message(glue('--- cell groups: {groups}'))    
    }

    if (length(cpg_groups(db)) > 0){
        groups <- paste(cpg_groups(db), collapse = ', ')
        message(glue('--- interval groups: {groups}'))    
    }

}

#' Load plpdb database from disk
#' 
#' @param db_root root directory of the cgdb database
#' 
#' @export
plpdb_load <- function(db_root){
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
    db <- new('plpdb', db_root = db_root, cpgs = cpgs, cells=cells, CPG_NUM=max(cpgs$id))    
    
    return(db)    
}


#' Add plate to plpdb frpm smat reads
#' 
#' @param prefix path prefix of the location of the smat files
#' 
#' @export
plpdb_add_plate <- function(db, prefix, plate_name=NULL, overwrite=TRUE, update_cells=TRUE, verbose=TRUE){
    reads_dir <- paste0(prefix, '_reads')
    if (!dir.exists(reads_dir)){
        stop("reads dir doesn't exist")
    }

    cell_metadata_file <- paste0(prefix, '_cell_metadata.tsv')
    cell_metadata <- fread(cell_metadata_file)
    smat_colnames <- fread( paste0(prefix, '_colnames.tsv')) %>% as.tibble() %>% pull(1)

    cells <- cell_metadata %>% 
        left_join(tibble(cell_id = smat_colnames), by = "cell_id") %>% 
        select(-one_of('plate')) %>% 
        separate(cell_id, c('plate', 'cell_num'), sep='\\.', remove=FALSE) %>%
        arrange(cell_num) %>% 
        as_tibble()
    
    plyr::a_ply(cells, 1, function(x) plpbdb_add_cell(db=db, reads_dir=reads_dir, cell_id = x$cell_id, cell=x$cell_num, plate=x$plate, overwrite=overwrite), .parallel = TRUE)

    if (update_cells){
        db@cells <- fread(glue('{db@db_root}/cells.csv')) %>% as_tibble()
        db <- cgdb_update_cells(db, cells, append=TRUE)    
    } else {
        warning('cells.csv not updated')
    }    
    return(db)
}

plpbdb_add_cell <- function(db, reads_dir, cell_id, cell, plate, overwrite){
    reads <- fread(glue("{reads_dir}/{cell_id}.tsv")) %>% select(chrom, start) %>% filter(chrom %in% gintervals.all()$chrom) %>% mutate(end = start + 1) %>% gpatterns:::.gpatterns.force_chromosomes() %>% as_tibble()

    pileup <- pileup_cell_reads(db, reads) %>% mutate(meth = 0)    
    pileup <- pileup %>% mutate(cov = as.numeric(cov))
    db <- cgdb_add_cell(db, pileup, cell, plate, overwrite=overwrite)
    return(db)
}

pileup_cell_reads <- function(db, reads){
    opt <- options(gmax.data.size=1e9)
    on.exit(options(opt))
    pileup <- reads %>% gintervals.neighbors1(db@cpgs) %>% filter(dist == 0) %>% group_by(chrom1, start1, end1, id) %>% summarise(cov = n()) %>% ungroup() %>% rename(chrom=chrom1, start=start1, end=end1) %>% arrange(id)

    return(pileup)
}

cells_pileup <- function(db, binsize){
    db %>% summarise_intervals(giterator.intervals(iterator=binsize)) %>% select(-meth)
}

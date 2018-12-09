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
        df <- df %>% inner_join(db@cpgs %>% select(chrom, start, end, id) ,by=c('chrom', 'start', 'end')) %>% filter(!is.na(id)) %>% arrange(id)

        cov_vec <- as.numeric(df$cov)
        met_vec <- as.numeric(df$meth)
        idxs <- as.integer(df$id)

        writeBin(idxs, paste0(fname, '.idx.bin'), size=4)                
        writeBin(met_vec, paste0(fname, '.meth.bin'), size=4)
        writeBin(cov_vec, paste0(fname, '.cov.bin'), size=4)        
        message(glue('created {cell}'))    
    }
    invisible(db)
}

lock_db <- function(db){

    lock <- filelock::lock(glue('{db@db_root}/.cells_lock'), timeout=Inf)
    return(lock)
}

unlock_db <- function(lock){

    filelock::unlock(lock)
}

cgdb_update_cells <- function(db, cells, append=FALSE){   
    l <- lock_db(db)
    db@cells <- fread(glue('{db@db_root}/cells.csv')) %>% as_tibble()

    if (append){
        if (!is.character(db@cells$cell_num)){
            db@cells$cell_num <- as.character(db@cells$cell_num)
        }
        if (!is.character(cells$cell_num)){
            cells$cell_num <- as.character(cells$cell_num)
        }
        cells <- bind_rows(db@cells %>% filter(!(cell_id %in% cells$cell_id)), cells)
    }
    
    fwrite(cells, glue('{db@db_root}/cells.csv'), sep=',')    
    unlock_db(l)
    db@cells <- cells
    invisible(db)
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
        arrange(cell_num)

    if (has_stats(smat)){        
        cells <- cells %>% left_join(smat$stats)
    }

    message('Converting smat to data frame')
    smat_df <- smat.to_df(smat) %>% select(-id)    
    db <- cgdb_add_plate_from_df(db, smat_df, cells, plate_name=plate_name, overwrite=overwrite, update_cells=update_cells, verbose=verbose)
   
    return(db)
}


#' Add data from tidy_cpgs directory to cgdb
#' 
#' @param cell_id string with the cell id, needs to be of the format "{plate}.{num}", where num has to be numeric
#' @param tidy_cpgs_dir path to the directory of the tidy_cpgs (unique)
#' @param db cgdb object
#' @param tidy_cpgs_dir_non_unique path to the directory of the non-unique tidy_cpgs for statistics computation (optionsl)
#' @param metadata data frame with cell_id and additional fields with cell metadata (optional)
#' @param single_cell take only a single methylation call per CpG (chosen randomly)
#' @param overwrite overwrite existing cells
#' @param update_cells update cells metadata for the db object
#' 
#' @export
cgdb_add_cell_from_tidy_cpgs_dir <- function(cell_id, tidy_cpgs_dir, db, tidy_cpgs_dir_non_unique=NULL, metadata=NULL, single_cell=TRUE, overwrite=TRUE, update_cells=TRUE){
    tidy_cpgs <- gpatterns:::.gpatterns.get_tidy_cpgs_from_dir(tidy_cpgs_dir, uniq=TRUE)

    if (!is.null(tidy_cpgs_dir_non_unique)){
        tidy_cpgs_stats_dir <- paste0(tidy_cpgs_dir_non_unique, '/stats')
        uniq_tidy_cpgs_stats_dir <- paste0(tidy_cpgs_dir, '/stats')
        stats <- gpatterns::gpatterns.get_tcpgs_stats(tidy_cpgs_stats_dir, uniq_tidy_cpgs_stats_dir)$stats 
        metadata <- metadata %>% bind_cols(stats)
    } 

    cgdb_add_cell_from_tidy_cpgs(cell_id=cell_id, tidy_cpgs=tidy_cpgs, metadata=metadata, db=db, single_cell=single_cell, overwrite=overwrite, update_cells=update_cells)

}

#' Add data from tidy_cpgs data frame to cgdb
#' 
#' @param tidy_cpgs tidy cpgs data frame
#' 
#' @inheritParams cgdb_add_cell_from_tidy_cpgs_dir
#' 
#' @export
cgdb_add_cell_from_tidy_cpgs <- function(cell_id, tidy_cpgs, db, metadata=NULL, single_cell=TRUE, overwrite=TRUE, update_cells=TRUE){

    if (single_cell){
        calls <- tcpgs2calls(tidy_cpgs, cell_id) %>% 
            mutate(unmeth = if_else(meth == 0, 1, 0), cov=1) %>% 
            select(chrom, start, end, meth, unmeth, cov)
    } else {
        calls <- gpatterns::gpatterns.tidy_cpgs_2_pileup(tidy_cpgs) %>% mutate(cell_id=!!cell_id) %>% select(chrom, start, end, meth, unmeth, cov)
    }       
    
    cgdb_add_cell_from_df(df=calls, cell_id=cell_id, metadata=metadata, db=db, overwrite=overwrite, update_cells=update_cells)
}

parse_cell_id <- function(cell_id){
    parts <- stringr::str_split(cell_id, '\\.')[[1]]
    if (length(parts) < 2){
        stop('cell_id needs to be in the form of {plate}.{num}')
    }

    cell_num <- as.numeric(parts[length(parts)] )
    if (is.na(cell_num)){
        stop('cell_id needs to be in the form of {plate}.{num}')   
    }

    plate <- paste(parts[-length(parts)], collapse='.')
    return(list(plate=plate, cell_num=cell_num))
}

#' Add data from data frame to cgdb
#' 
#' @param df data frame with the following fields: 
#' \itemize{
#'  \item{"cell_id" }{string with the cell id, needs to be of the format "{plate}.{num}", where num has to be numeric}
#'  \item{"chrom" }{chromosome, (character)}
#'  \item{"start" }{start position of the CpG (numeric)}
#'  \item{"end" }{end position of the CpG (start + 1, numeric)}
#'  \item{"cov" }{number of times the CpG was covered (numeric)}
#'  \item{"meth" }{number of times the CpG was methylated (numeric)}
#'  \item{"unmeth" }{number of times the CpG was unmethylated (cov - meth, numeric)}
#' }
#' 
#' @inheritParams cgdb_add_cell_from_tidy_cpgs_dir
#' 
#' @export
cgdb_add_cell_from_df <- function(df, db, cell_id, metadata=NULL, overwrite=TRUE, update_cells=TRUE){
    parsed_cell_id <- parse_cell_id(cell_id)
    plate <- parsed_cell_id$plate
    cell_num <- parsed_cell_id$cell_num    
    cgdb_add_cell(db=db, df=df, cell=cell_num, plate=plate, overwrite=overwrite)

    if (update_cells){
        metadata <- bind_cols(tibble(cell_id=cell_id, plate=plate, cell_num=cell_num), metadata)        
        db <- cgdb_update_cells(db, metadata, append=TRUE)    
    } else {
        warning('cells.csv not updated')
    }    

    invisible(db)        
}
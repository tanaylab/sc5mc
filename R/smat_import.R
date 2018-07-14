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
#' @param name name of the object
#' @param description description of the object
#' @param keep_tidy_cpgs keep tidy cpgs attached to the smat object
#' @param load_existing load existing tidy cpgs
#' @param filter_dups relevenat only when load_existing is TRUE: filter duplicates
#' @param single_cell take only a single methylation call per CpG
#' @param ... additional parameters to \code{gpatterns::gpatterns.import_from_bam}
#' 
#' @return smat object (inivisibly if prefix is not NULL)
#'
#' @export
smat.from_bams <- function(metadata, groot, prefix=NULL, workdir=tempdir(), use_sge = TRUE, name='', description='', keep_tidy_cpgs=TRUE, load_existing=FALSE, filter_dups=FALSE, io_saturation=FALSE, threads=10, cell_metadata=NULL, keep_tidy_calls=TRUE, single_cell=TRUE, ...){

    bam2calls <- function(bams, lib, workdir=workdir, keep_tidy_cpgs=keep_tidy_cpgs, load_existing=load_existing, ...){                        
        if (keep_tidy_cpgs){
            track_workdir <- paste0(prefix, '_tcpgs/', lib)            
        } else {
            track_workdir <- tempfile(tmpdir=workdir)            
            on.exit(system(sprintf('rm -rf %s', track_workdir)))            
        }
        system(sprintf('mkdir -p %s', track_workdir))
        
        if (!load_existing){            
            gpatterns::gpatterns.import_from_bam(bams, workdir=track_workdir, steps=c('bam2tidy_cpgs', 'filter_dups'), groot=groot, ...)    
        } else {
            if (filter_dups){
                gpatterns::gpatterns.import_from_bam('', workdir=track_workdir, steps=c('filter_dups'), groot=groot, ...)    
            }
        }

        tcpgs_dir <- paste0(track_workdir, '/tidy_cpgs')
        uniq_tcpgs_dir <- paste0(track_workdir, '/tidy_cpgs_uniq')
        tcpgs <- gpatterns:::.gpatterns.get_tidy_cpgs_from_dir(uniq_tcpgs_dir)        

        if (is.null(tcpgs)){
            return(NULL)
        }
        if (single_cell){
            calls <- tcpgs2calls(tcpgs, lib)    
        } else {
            calls <- gpatterns::gpatterns.tidy_cpgs_2_pileup(tcpgs) %>% mutate(track=lib) %>% select(track, chrom, start, end, meth, unmeth, cov)
        }        

        tidy_cpgs_stats_dir <- paste0(tcpgs_dir, '/stats')
        uniq_tidy_cpgs_stats_dir <- paste0(uniq_tcpgs_dir, '/stats')
        stats <- gpatterns::gpatterns.get_tcpgs_stats(tidy_cpgs_stats_dir, uniq_tidy_cpgs_stats_dir) %>% .$stats %>% mutate(lib = lib) %>% select(lib, everything())
        stats[['cpg_num']] <- nrow(calls)       

        return(list(calls=calls, stats=stats))
    }

    gsetroot(groot)    

    cmds <- plyr::daply(metadata, plyr::.(lib), function(x) {
        bams <- paste(qqv("'@{x$bam}'"), collapse=', ')
        qq("bam2calls(c(@{bams}), lib = '@{x$lib[1]}', workdir='@{workdir}', keep_tidy_cpgs=@{keep_tidy_cpgs}, load_existing=@{load_existing}, ...)")
    })
    
    
    if (use_sge){        
        res <- gpatterns::gcluster.run2(command_list=cmds, io_saturation=io_saturation, threads=threads, packages='gpatterns')
        tidy_calls <- res %>% map('retv') %>% map('calls') %>% compact() %>% map_df(~ .x)
        stats <- res %>% map('retv') %>% map('stats') %>% compact() %>% map_df(~ .x)
    } else {
        tidy_calls <- map(cmds, ~ eval(parse(text=.x))) %>% map('calls') %>% compact() %>% map_df(~ .x)
        stats <- map(cmds, ~ eval(parse(text=.x))) %>% map('stats') %>% compact() %>% map_df(~ .x)
    }
    
    if (single_cell){       
        if (nrow(tidy_calls) == 0) {
            logerror('did not create any tidy calls')
        }
        tidy_calls <- tidy_calls %>% mutate(unmeth = if_else(meth == 0, 1, 0), cov=1)
    }
    if (keep_tidy_calls){
        fwrite(tidy_calls, paste0(prefix, '_tidy_calls.csv'))    
    }    
    
    smat <- .tidy_calls2smat(tidy_calls)
    smat$stats <- stats
    smat$name <- name
    smat$description <- description

    if (!is.null(smat$stats)){
        smat$stats <- smat.cell_marginals(smat) %>% rename(lib = cell_id, cg_num = cov) %>% inner_join(smat$stats)    
    }
    

    if (!is.null(cell_metadata)){
        smat$cell_metadata <- fread(cell_metadata[1])
        smat$stats <- smat$cell_metadata %>% full_join(smat$stats %>% rename(cell_id = lib))
    }

    if (keep_tidy_cpgs){
        smat$tidy_cpgs <- paste0(prefix, '_tcpgs') 
    }

    if (!is.null(prefix)){
        system(qq('mkdir -p @{dirname(prefix)}'))
        smat.save(smat, prefix)
        invisible(smat)
    }
    return(smat)
}

# #' @export
# smat.from_tidy_calls <- function(tidy_calls, ){
#     smat.from_bams(data.frame(lib = cells, bam=''), prefix=prefix, keep_tidy_cpgs=TRUE, load_existing=TRUE, ...)
# }


#' Create smat object from misha tracks
#'
#' @param tracks names of tracks
#'
#' @param libs names of new object smat columns
#' @param prefix path prefix in which to save the smat files. if NULL the object
#' would not be written to disk
#' @param use_sge use sun grid engine cluster
#' @param name name of the object
#' @param description description of the object
#' @param ... additional parameters to gcluster.run2
#'
#' @return smat object
#'
#' @export
smat.from_tracks <- function(tracks, libs, prefix=NULL, use_sge=TRUE, name='', description='', ...){
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

    smat$name <- name
    smat$description <- description

    if (!is.null(prefix)){
        smat.save(smat, prefix)
    }

    return(smat)
}

#' Create smat object from data frame
#' @param df intervals ('chrom', 'start', 'end' fields) with additional fields:
#' meth (methylated calls), unmeth (unmethylated calls) and cov (total coverage)
#' @param name name of the object
#' @param description description of the object
#' @param intervs intervals set to take
#' 
#' @return smat object
#' @export
smat.from_df <- function(df, name='', description='', intervs=NULL){
    # .tidy_calls2smat(df %>% unite('coord', chrom, start, end), column_name='cell')
    smat <- .tidy_calls2smat(df, column_name='cell', intervs=intervs)
    smat$name <- name
    smat$description <- description    
    return(smat)
}


# Utils
.tidy_calls2smat <- function(tidy_calls, column_name='track', intervs=NULL){             
    message('creating intervs')  
    if (!is.null(intervs)){
        intervs <- intervs %>% arrange(chrom, start, end) %>% mutate(id = 1:n())
    } else {
        intervs <- tidy_calls %>% arrange(chrom, start, end) %>% distinct(chrom, start, end) %>% mutate(id = 1:n())    
    }
    
    
    if (has_name(tidy_calls, 'id')){
        tidy_calls <- tidy_calls %>% select(-id)
    }
    tidy_calls <- tidy_calls %>% inner_join(intervs, by=c('chrom', 'start', 'end')) %>% arrange(id)

    # parallel <- TRUE
    # if (nrow(tidy_calls) * ncol(tidy_calls) > (2^31 - 1)){
        parallel <- FALSE
    # }
    smat <- plyr::alply(c('meth', 'unmeth', 'cov'), 1, function(x) {
        message(sprintf('creating %s', x))
        tidy2smat(tidy_calls, 'id', column_name, x)
    }, .parallel=FALSE)

    names(smat) <- c('meth', 'unmeth', 'cov')
      
    smat$intervs <- intervs

    for (.x in c('meth', 'unmeth', 'cov')){        
        stopifnot(all(rownames(smat[[.x]]) == smat$intervs$id))
        # rownames(smat[[.x]]) <- smat$intervs$id
    }

    class(smat) <- 'smat'

    return(smat)
}

#' Create sparse matrix from tidy data frame
tidy2smat <- function(data, row, column, value, ...){
    row_u <- unique(data[[row]])
    i <- match(data[[row]], row_u)

    col_u <- unique(data[[column]])
    j <- match(data[[column]], col_u)

    val <- data[[value]]

    Matrix::sparseMatrix(i = i, j = j, x = val, dimnames = list(as.character(row_u), as.character(col_u)), ...)
}

tcpgs2calls <- function(tcpgs, track){
   tcpgs %>%
        arrange(chrom, cg_pos) %>%
        distinct(chrom, cg_pos, .keep_all=T) %>%
        mutate(start=cg_pos, end=start+1, track=track) %>%
        select(track, chrom, start, end, meth)
}


# smat.merge_by_reads <- function(smat_list, do_merge=TRUE, filter_dups=TRUE, prefix=NULL, workdir=tempdir(), groot=GROOT, ...){
#     if (!is.null(prefix)){
#          workdir <- paste0(prefix, '_tcpgs/')
#          system(glue('mkdir -p {workdir}'))
#     }
#     metadata <- map(smat_list, get_smat_tcpgs_files)
#     metadata <- map_dfr(metadata, ~ .x %>% tidyr::unnest(files) ) %>% mutate(fn = basename(files)) %>% group_by(cell_id, fn) %>% nest(files) %>% mutate(out_dir = glue('{workdir}/{cell_id}/tidy_cpgs')) %>% mutate(id = 1:n())
    
#     merge_tcpgs <- function(i){
#         # print(scales::percent(i / nrow(metadata)))
#         .x <- metadata[i, ]
#         system(glue('mkdir -p {.x$out_dir}'))
#         files <- .x$data[[1]][[1]]
#         system(glue('cp {files[1]} {.x$out_dir}/{.x$fn}'))
#         walk(files[-1], function(f) system(glue('gzip -d -c {f} | tail -n +2 | gzip -c >> {.x$out_dir}/{.x$fn}')))
#     }

#     if (do_merge){
#         message('merging tidy cpgs')    
#         res <- dlply(metadata, .(cell_id), function(x) {print(x$cell_id[1]); walk(x$id, merge_tcpgs)}, .parallel=TRUE)    
#     }
#     # a <- alply(1:nrow(metadata), 1, merge_tcpgs, .parallel=TRUE)
#     metadata <- metadata %>% distinct(cell_id) %>% rename(lib=cell_id) %>% mutate(bam = '')

#     #TODO: merge cell metadata and stats!
#     smat <- smat.from_bams(metadata, groot, prefix=prefix, workdir=workdir, load_existing=TRUE, filter_dups=TRUE, ...)
#     browser()
# }


# get_smat_tcpgs_files <- function(smat){
#     tibble(cell_id = list.files(smat$tidy_cpgs_dir)) %>% 
#         purrrlyr::by_row(~ list.files(glue('{smat$tidy_cpgs_dir}/{.x$cell_id}/tidy_cpgs'), pattern='.*\\.tcpgs.gz$', full.names=TRUE), .to='files')
# }

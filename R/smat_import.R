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
#' @return smat object (inivisibly if prefix is not NULL)
#'
#' @export
smat.from_bams <- function(metadata, groot, prefix=NULL, workdir=tempdir(), use_sge = TRUE, ...){
    bam2calls <- function(bams, lib, workdir=workdir, use_sge=TRUE){                  
        track_workdir <- tempfile(tmpdir=workdir)
        system(sprintf('mkdir -p %s', track_workdir))
        on.exit(system(sprintf('rm -rf %s', track_workdir)))        
        gpatterns::gpatterns.import_from_bam(bams, workdir=track_workdir, steps=c('bam2tidy_cpgs', 'filter_dups'), groot=groot)

        tcpgs_dir <- paste0(track_workdir, '/tidy_cpgs')
        uniq_tcpgs_dir <- paste0(track_workdir, '/tidy_cpgs_uniq')
        tcpgs <- .gpatterns.get_tidy_cpgs_from_dir(uniq_tcpgs_dir)        

        if (is.null(tcpgs)){
            return(NULL)
        }
        calls <- tcpgs2calls(tcpgs, lib)

        tidy_cpgs_stats_dir <- paste0(tcpgs_dir, '/stats')
        uniq_tidy_cpgs_stats_dir <- paste0(uniq_tcpgs_dir, '/stats')
        stats <- gpatterns.get_tcpgs_stats(tidy_cpgs_stats_dir, uniq_tidy_cpgs_stats_dir) %>% .$stats %>% mutate(lib = lib) %>% select(lib, everything())
        stats[['cpg_num']] <- nrow(cols)       

        return(list(calls=calls, stats=stats))
    }

    gsetroot(groot)    

    cmds <- plyr::daply(metadata, .(lib), function(x) {
        bams <- paste(qqv("'@{x$bam}'"), collapse=', ')
        qq("bam2calls(c(@{bams}), lib = '@{x$lib[1]}', workdir='@{workdir}')")
    })
    
    if (use_sge){        
        res <- gcluster.run2(command_list=cmds, io_saturation=T, packages='gpatterns')
        tidy_calls <- res %>% map('retv') %>% map('calls') %>% compact() %>% map_df(~ .x)
        stats <- res %>% map('retv') %>% map('stats') %>% compact() %>% map_df(~ .x)
    } else {
        tidy_calls <- map(cmds, ~ eval(parse(text=.x))) %>% map('calls') %>% compact() %>% map_df(~ .x)
        stats <- map(cmds, ~ eval(parse(text=.x))) %>% map('stats') %>% compact() %>% map_df(~ .x)
    }

    tidy_calls <- tidy_calls %>% mutate(unmeth = if_else(meth == 0, 1, 0), cov=1)

    smat <- .tidy_calls2smat(tidy_calls)
    smat$stats <- stats

    if (!is.null(prefix)){
        system(qq('mkdir -p @{dirname(prefix)}'))
        smat.save(smat, prefix)
        invisible(smat)
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
    .tidy_calls2smat(df %>% unite('coord', chrom, start, end), column_name='cell')
}


# Utils

.tidy_calls2smat <- function(tidy_calls, column_name='track'){
    smat <- plyr::alply(c('meth', 'unmeth', 'cov'), 1, function(x) {
        message(sprintf('creating %s', x))
        tidy2smat(tidy_calls, 'coord', column_name, x)
    }, .parallel=TRUE)

    names(smat) <- c('meth', 'unmeth', 'cov')

    message('creating intervs')
    smat$intervs <- tibble(coord = rownames(smat[[1]])) %>% 
        separate(coord, c('chrom', 'start', 'end')) %>% mutate(id = 1:n()) %>% 
        mutate(start = as.numeric(start), end = as.numeric(end))

    return(smat)
}

#' Create sparse matrix from tidy data frame
tidy2smat <- function(data, row, column, value, ...){
    row_u <- unique(data[[row]])
    i <- match(data[[row]], row_u)

    col_u <- unique(data[[column]])
    j <- match(data[[column]], col_u)

    val <- data[[value]]    

    Matrix::sparseMatrix(i = i, j = j, x = val, dimnames = list(row_u, col_u), ...)
}

tcpgs2calls <- function(tcpgs, track){
   tcpgs %>%
        arrange(chrom, cg_pos) %>%
        distinct(chrom, cg_pos, .keep_all=T) %>%
        mutate(start=cg_pos, end=start+1, track=track) %>%
        select(track, chrom, start, end, meth) %>%
        unite('coord', chrom:end)
}
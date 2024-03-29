#' @export
downsample_tab <- function(tab, n = NULL){

    if (!is.null(rownames(tab))){
        rownames_tab <- rownames(tab)
    } else {
        rownames_tab <- NULL
    }

    tab <- as_tibble(tab)
    if (is.null(n)){
        n <- min(colSums(tab))
    }
    for (i in 1:ncol(tab)){
        v <- tab[[i]]
        tab[[i]] <- tabulate(sample(rep(1:length(v),times=v), replace=FALSE,size=n), nbins=length(v))
    }

    if (!is.null(rownames_tab)){
        tab <- as.data.frame(tab)
        rownames(tab) <- rownames_tab
    }

    return(tab)
}

#' @export
downsample_df <- function(df, columns, n=NULL){    
    tab <- df %>% select(!!! rlang::syms(columns))        
    ds_tab <- downsample_tab(t(tab), n=n)
    ds_tab <- t(ds_tab) %>% as.data.frame() %>% rlang::set_names(colnames(tab))
    return(df %>% select(-one_of(columns)) %>% bind_cols(ds_tab))
}

#' @export
downsample_cpgs <- function(tab, n=NULL){
	tab <- tab %>% select(cell_id, cov, meth)
	if (is.null(n)){
        n <- min(colSums(tab))
    }

    ds_tab <- tab %>% 
    	filter(cov >= n) %>% 
    	mutate(unmeth = cov - meth) %>% 
    	select(cell_id, meth, unmeth) %>% 
    	as.data.frame() %>% 
    	tibble::column_to_rownames('cell_id') %>% 
    	t() %>% 
    	downsample_tab(n=n) %>%
    	t() %>% 
    	as.data.frame() %>% 
    	tibble::rownames_to_column('cell_id') %>% 
    	rlang::set_names(c('cell_id', 'meth', 'unmeth')) %>% 
    	as_tibble()

    ds_tab <- ds_tab %>% mutate(cov = meth + unmeth) %>% select(cell_id, cov, meth)

    return(ds_tab)
}

#' @export
downsample_meth <- function(df, dsn=NULL, group_vars='cell_id', .parallel=TRUE){
    grouping <- rlang::syms(group_vars)
    if (is.null(dsn)){
        dsn <- df %>% group_by(!!!grouping) %>% summarise(n = sum(cov)) %>% pull(n) %>% min()
    }

    df <- df %>% select(cell_id, chrom, start, end, cov, meth) %>% mutate(unmeth = cov - meth)

    df_ds <- plyr::ddply(df, group_vars, function(x){       

        d_meth <- x %>% select(chrom, start, end, meth) %>% mutate(m = 1) %>% uncount(meth)
        d_unmeth <- x %>% select(chrom, start, end, unmeth) %>% mutate(m = 0) %>% uncount(unmeth)
            
        res <- bind_rows(d_meth, d_unmeth) %>% 
                sample_n(dsn) %>% 
                group_by(chrom, start, end) %>% 
                summarise(cov = n(), meth = sum(m)) %>% 
                ungroup()

        return(as.data.frame(res))
    }, .parallel = .parallel)

    df_ds <- as_tibble(df_ds)

    return(df_ds)
}


# debug utils
read_cell_binary_files <- function(db, plate, cell, n=1e6){    
    tibble(id = readBin(glue('{db@db_root}/data/{plate}/{cell}.idx.bin'), what=integer(), size=4, n=n), cov = readBin(glue('{db@db_root}/data/{plate}/{cell}.cov.bin'), what=numeric(), size=4, n=n), meth = readBin(glue('{db@db_root}/data/{plate}/{cell}.meth.bin'), what=numeric(), size=4, n=n))
}
downsample_tab <- function(tab, n = NULL){
    tab <- as_tibble(tab)
    if (is.null(n)){
        n <- min(colSums(tab))
    }
    for (i in 1:ncol(tab)){
        v <- tab[[i]]
        tab[[i]] <- tabulate(sample(rep(1:length(v),times=v), replace=FALSE,size=n), nbins=length(v))
    }
    return(tab)
}

downsample_df <- function(df, columns, n=NULL){    
    tab <- df %>% select(!!! rlang::syms(columns))        
    ds_tab <- downsample_tab(t(tab), n=n)
    ds_tab <- t(ds_tab) %>% as.data.frame() %>% rlang::set_names(colnames(tab))
    return(df %>% select(-one_of(columns)) %>% bind_cols(ds_tab))
}

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
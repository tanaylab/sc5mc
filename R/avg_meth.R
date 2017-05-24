#' @export
sc5mc.global_meth_trend <- function(smat, breaks=seq(0, 0.1, by=0.002), include.lowest=TRUE, min_cpgs=500, parallel=TRUE){
    gmeth <- smat.summarise_by_track(smat, gpatterns:::.gpatterns.cg_cont_500_track, breaks=breaks, include.lowest=include.lowest, group_name='cg_cont', parallel=parallel)
    breaks_numeric_df <- tibble(breaks_numeric=zoo::rollmean(breaks, k=2), cg_cont=levels(gmeth$cg_cont))
    gmeth <- gmeth %>%
    	mutate(cg_cont = as.character(cg_cont)) %>% 
    	left_join(breaks_numeric_df, by='cg_cont') %>% 
    	select(cg_cont, breaks_numeric, track, ncpgs, meth)

    gmeth <- gmeth %>% filter(ncpgs >= min_cpgs)
    return(gmeth)
}
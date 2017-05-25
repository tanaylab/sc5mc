#' @param smat smat object
#'
#' @param breaks breaks of cpg content
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param min_cpgs minimal number of CpGs per cg content strat
#' @param parallel compute parallely per group (using doMC package)
#'
#' @return data frame with the following fields:
#' cg_cont` with the CpG content bin,
#' `breaks_numeric` with the middle point between each bin
#' `cell` cell name
#' `ncpgs` number of CpGs in the bin
#' `meth` number of methylated CpGs
#' `unmeth` number of unmethylated CpGs
#' `avg` fraction of methylated calls out of all calls
#'
#'
#' @export
sc5mc.global_meth_trend <- function(smat,
                                    breaks = seq(0, 0.1, by = 0.002),
                                    include.lowest = TRUE,
                                    min_cpgs = 500,
                                    parallel = TRUE) {

    gmeth <- smat.summarise_by_track(
        smat,
        gpatterns:::.gpatterns.cg_cont_500_track,
        breaks = breaks,
        include.lowest = include.lowest,
        group_name = 'cg_cont',
        parallel = parallel
    )
    breaks_numeric_df <- tibble(breaks_numeric=zoo::rollmean(breaks, k=2), cg_cont=levels(gmeth$cg_cont))
    gmeth <- gmeth %>%
    	mutate(cg_cont = as.character(cg_cont)) %>%
    	left_join(breaks_numeric_df, by='cg_cont') %>%
    	select(cg_cont, breaks_numeric, cell, ncpgs, meth, unmeth, avg)

    gmeth <- gmeth %>% filter(ncpgs >= min_cpgs)
    return(gmeth)
}

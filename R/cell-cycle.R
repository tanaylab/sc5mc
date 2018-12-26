get_cell_tor_cov <- function(db, tor_breaks, tor_labels, tor_track=NULL){
	if (!has_name(db@cpgs, 'tor')){
		if (is.null(tor_track)){
			stop('Please provide tor_track')
		}
		gvtrack.create("tor", tor_track, "avg")
    	gvtrack.iterator("tor", sshif=-15000, eshift=15000)
    	opt <- options(gmax.data.size=1e9)
    	on.exit(options(opt))
    	db <- db %>% g_mutate_cpgs('tor')
	}
	db <- db %>% mutate_cpgs(repli = cut(tor, breaks=tor_breaks, labels=tor_labels))
	cell_tor <- db %>% group_by_cpgs(repli) %>% summarise() 

	return(cell_tor)

}


#' Phase cells (to G and S) using ratio of coverage in late and early replicating regions
#' 
#' @param db cgdb object
#' @param early_tor early_tor
#' @param late_tor late_tor
#' @param ratio_thresh ratio_thresh
#' @param tor_track tor_track
#' 
#' @examples
#' \dontrun{
#' 		 cell_phases <- phase_cells(db, late_tor=c(0, 50), early_tor=c(66, 83), ratio_thresh=1.42, tor_track='encode.repliseq.wgEncodeUwRepliSeqImr90WaveSignalRep1')
#' }
#' 
#'   
#' @export
phase_cells <- function(db, early_tor, late_tor, ratio_thresh, tor_track=NULL){
	cell_tor <- get_cell_tor_cov(db=db, tor_breaks=c(late_tor, early_tor), tor_labels=c('late_cov', 'mid_cov', 'early_cov'))
	cell_tor <- cell_tor %>% select(cell_id, repli, cov) %>% spread(repli, cov)   
	cell_tor <- cell_tor  %>% mutate(phase = ifelse((early_cov / late_cov) >= ratio_thresh, 'S', 'G')) 	

	return(cell_tor)
}
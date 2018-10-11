#' Metacell layout using force directed projection of a low degree mc graph
#'
#' @param mc2d_id 2d object to add
#' @param mc_id meta cell id to work with
#' @param graph_id graph_id of the similarity graph on cells from the metacell
#' @param mc_subset a subset of metacells to project (NULL by default)
#'
#' @export
sc5mc.mc2d_force_knn <- function(mc_map, cgraph)
{
	mc <- mc_map$mc  
	mc <- as.numeric(as.factor(mc))
    names(mc) <- mc_map$cell_id
	
	mgraph <- mc5mc_comp_mgraph(mc, cgraph)	
	cl_xy <- mc5mc_comp_graph_coord(mgraph, N=max(mc, na.rm=T))
	xy <- mc5mc_comp_cell_coord(mc, cgraph, mgraph, cl_xy)	
	mc_names <- levels(factor(mc_map$mc))

	cell_coord <- as_tibble(xy) %>% 
		mutate(cell_id = names(xy$x)) %>% 
		select(cell_id, x, y) %>% 
		left_join(mc_map)

	mc_coord <- as_tibble(cl_xy[c('x_cl', 'y_cl')]) %>% 
		mutate(mc = mc_names[as.numeric(names(cl_xy$x_cl))]) %>% 
		select(mc, x=x_cl, y=y_cl)

	edges <- mgraph %>% 
		mutate(mc1 = mc_names[mc1], mc2 = mc_names[mc2]) %>% 
		left_join(mc_coord %>% rename(mc1 = mc, x_mc1 = x, y_mc1 = y), by='mc1')  %>% 
		left_join(mc_coord %>% rename(mc2 = mc, x_mc2 = x, y_mc2 = y), by='mc2') %>% 
		filter(mc1 != mc2) %>%
		as_tibble()

	return(list(cells=cell_coord, mc=mc_coord, edges=edges, graph=cl_xy$g))
}

mc_confusion_mat <- function(mc, graph)
{
	mc_map = rep(NA, length(levels(graph$cell1)))
	mc_map[factor(names(mc),levels(graph$cell1))] = mc
	mc1 = mc_map[graph$cell1]
	mc2 = mc_map[graph$cell2]
	N = max(mc1, mc2, na.rm=T)

	confu = matrix(tabulate((mc1-1)*N+mc2, nbins=N*N), nrow=N)
	rownames(confu) = 1:N
	colnames(confu) = 1:N
	return(confu)
}


#' @export
mc5mc_comp_mgraph <- function(mc, graph)
{
	mc2d_K = get_param("mcell_mc2d_K")
	mc2d_T_edge = get_param("mcell_mc2d_T_edge")
	mc2d_max_confu_deg = get_param("mcell_mc2d_max_confu_deg")
	mc2d_edge_asym = get_param("mcell_mc2d_edge_asym")

	restrict_in_degree = TRUE

	confu = mc_confusion_mat(mc, graph)

	csize = as.matrix(table(mc))
	csize = pmax(csize, 20)
	csize2 = csize %*% t(csize)
	csize2 = csize2 / median(csize)**2
	confu = confu / csize2

	confu_p_from = confu/rowSums(confu)
	confu_p_to = t(confu)/colSums(confu)

	if(!is.null(mc2d_max_confu_deg)) {
		rank_fr = t(apply(confu_p_from, 1, rank))
		rank_to = t(apply(confu_p_to, 1, rank))
		rank2 = rank_fr * rank_to
		diag(rank2) = 1e+6
		amgraph = apply(rank2, 1, function(x) {  rank(-x) <= (1+mc2d_max_confu_deg) })
		mgraph = amgraph * ((confu_p_from + confu_p_to) > mc2d_T_edge)
		if(restrict_in_degree) {
			amgraph2 = t(apply(rank2, 2, function(x) {  rank(-x) <= (1+mc2d_max_confu_deg) }))
			mgraph = mgraph * amgraph2
		}

		if(mc2d_edge_asym) {
			mgraph = amgraph * (confu_p_from>mc2d_T_edge)
			mgraph = amgraph * (t(confu_p_to)>mc2d_T_edge)
		}
		mgraph = mgraph>0 | t(mgraph>0)
	} else {
		mgraph = (confu_p_from + confu_p_to) > mc2d_T_edge
	}

	N = nrow(mgraph)
	e = which(mgraph>0)
	n1 = ceiling((e)/N)
	n2 = 1+((e-1) %% N)
	return(data.frame(mc1 = n1, mc2 = n2))
}

#' @export
mc5mc_comp_graph_coord <- function(mc_graph, N)
{
	n1 = mc_graph$mc1
	n2 = mc_graph$mc2
	rEG <- new("graphNEL", nodes=as.character(1:N), edgemode="undirected")

	rEG = addEdge(as.character(n1[n1!=n2]), as.character(n2[n1!=n2]), rEG, rep(1, length(n1[n1!=n2])))

	g = layoutGraph(rEG, layoutType="neato")
	x_cl = nodeRenderInfo(g)$nodeX
	y_cl = nodeRenderInfo(g)$nodeY
	names(x_cl) = 1:N
	names(y_cl) = 1:N
	return(list(g=g, x_cl=x_cl, y_cl=y_cl))
}

#' @export
mc5mc_comp_cell_coord <- function(mc, graph, mgraph, cl_xy)
{
	mc2d_proj_blur = get_param("mcell_mc2d_proj_blur")
	mc2d_K_cellproj = get_param("mcell_mc2d_K_cellproj")

	x_cl = cl_xy$x_cl
	y_cl = cl_xy$y_cl

	N_mc = length(x_cl)
	N_c = length(mc)

	blurx = mc2d_proj_blur*(max(x_cl) - min(x_cl))
	blury = mc2d_proj_blur*(max(y_cl) - min(y_cl))

	is_active = rep(FALSE, N_mc*N_mc)
	is_active[(mgraph$mc1-1) * N_mc + mgraph$mc2] = TRUE
	is_active[((1:N_mc)-1) * N_mc + 1:N_mc] = TRUE

	mc_key1 = mc[levels(graph$cell1)]
	mc_key2 = mc[levels(graph$cell2)]
	mc1 = mc_key1[graph$cell1]
	mc2 = mc_key2[graph$cell2]

	f_in_mc = !is.na(mc1) & !is.na(mc2) #missing mc's, for example orphans
	f_active = is_active[(mc1-1)*N_mc + mc2]
	f = !is.na(f_active) & f_in_mc & f_active

#	deg = nrow(graph[f,])/length(mc)
#	T_w = 1-(mc2d_K_cellproj+1)/deg
#	f = f & graph$w > T_w

	to_x = x_cl[mc2]
	to_y = y_cl[mc2]
	c_x = tapply(to_x[f], graph$cell1[f], mean)
	c_y = tapply(to_y[f], graph$cell1[f], mean)

	base_x = min(c_x)
	base_y = min(c_y)
	max_x = max(c_x)
	base_x = base_x - (max_x-base_x)*0.1

	n_miss = length(setdiff(names(mc), names(c_x)))
	if(n_miss > 0) {
		stop("Missing coordinates in some cells that are not ourliers or ignored - check this out! (total ", n_miss, " cells are missing, maybe you used the wrong graph object?")
	}
	x = c_x[names(mc)]
	y = c_y[names(mc)]
	x = x + rnorm(mean=0, sd=blurx, n=length(x))
	y = y + rnorm(mean=0, sd=blury, n=length(y))

	return(list(x=x, y=y))
}

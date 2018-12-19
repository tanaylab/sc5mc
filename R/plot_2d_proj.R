#' Plot 2d projection of metacells
#' 
#' @param proj2d output of sc5mc.mc2d_force_knn
#' 
#' @export
sc5mc.plot_proj2d <- function(proj2d, cell_annot=NULL, mc_annot=NULL, mc_color_column=NULL, color_column=NULL, plot_mc_names = TRUE, mc_text_color='black', mc_text_size=2, mc_point_color='black', mc_point_fill='white', graph_color='black', graph_width = 0.5, point_size=0.5, point_color='black', alpha=1, mc_alpha=1, mc_point_size=5, plot_graph=TRUE, plot_mc_points=TRUE){
    if (!is.null(cell_annot)){
        proj2d$cells <- proj2d$cells %>% left_join(cell_annot)
        if (is.null(color_column)){
            stop('Please provide the name of the color column')
        }
        ggp <- proj2d$cells %>% ggplot(aes_string(x='x', y='y', color=color_column))
    } else {
        ggp <- proj2d$cells %>% ggplot(aes(x=x, y=y))
    }    

    if (is.null(color_column)){
        ggp <- ggp + 
            geom_point(size=point_size, alpha=alpha, color=point_color, shape=19)    
    } else {
        ggp <- ggp + 
            geom_point(size=point_size, alpha=alpha, shape=19)
    }
    

    if (plot_graph){
      ggp <- ggp + geom_segment(data=proj2d$edges, aes(x=x_mc1, y=y_mc1, xend=x_mc2, yend=y_mc2), color=graph_color, size=graph_width)   
    } 

    if (plot_mc_points){
        if (!is.null(mc_annot)){
            proj2d$mc <- proj2d$mc %>% left_join(mc_annot)
            if (is.null(mc_color_column)){
                stop('Please provide mc_color_column')
            }
            ggp <- ggp +  geom_point(data=proj2d$mc, color=mc_point_color, aes_string(fill=mc_color_column), size=mc_point_size, alpha=mc_alpha, shape=21)
        } else {
            ggp <- ggp + geom_point(data=proj2d$mc, color=mc_point_color, fill=mc_point_fill, size=mc_point_size, alpha=mc_alpha, shape=21)
        }    
    }


        
    if (plot_mc_names){
        ggp <- ggp + 
            geom_text(data = proj2d$mc, aes(label=mc), color=mc_text_color, size=mc_text_size)     
    }

    ggp <- ggp + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
   

    return(ggp)
}
#' plotYvalue function
#'
#' This function allows you to calculate and display the maximum and minimum values of expression.
#'
#' @param nr Number of components
#' @param max_color Maximum value of each component
#' @param min_color Minimum value of each component
#'
#' @keywords plotYvalue
#' @export
#'
plotYvalue <-
function(nr, max_color, min_color){

  pushViewport(viewport(layout = grid.layout(nr,1)))
  for(i in 1:nr){
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    grid.text(sprintf("%.1f", max_color[i]),
              x = 0.05, y= 0.95,
              hjust = 0,
              vjust = 1,
              gp = gpar(cex = 0.6)
    )
    grid.text(sprintf("%.1f", min_color[i]),
              x = 0.05, y= 0.05,
              hjust = 0,
              vjust = 0,
              gp = gpar(cex = 0.6)
    )
    popViewport()
  }
  popViewport()
}

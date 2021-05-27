#' plotColorKey2 function
#'
#' This function allows you to add a label about expressions of PAs.
#'
#' @param colors select color
#' @param min_color the minimum value of PA expression
#' @param max_color the maximum value of PA expression
#'
#' @keywords plotColorKey2
#' @export
#'
plotColorKey2 <-
function(colors, min_color, max_color){
  pushViewport(viewport(layout = grid.layout(
    1,2,
    widths = unit(c(4,1), c("null","null"))
  )))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))


  grid.raster(matrix(colors,nrow=1), y = 0.25, interpolate = F, width = 1, height = 0.5, vjust = 0)
  grid.rect(x=0, y = 0.25, width=1, height = 0.5, gp = gpar(fill = NA, col="black"), vjust = 0, hjust = 0)

  # grid.text(min_color, x=-0.1, y=0.24, hjust = 0, vjust = 1, gp = gpar(cex = 0.4))
  # if(max_color > 0 & min_color < 0){
  #   grid.text(0, x=(0-min_color)/(max_color-min_color), y=0.24, hjust = 0, vjust = 1, gp = gpar(cex = 0.4))
  # }
  grid.text(round(min_color,2), x=0, y=0.24, hjust = 0, vjust = 1, gp = gpar(cex = 0.8))

  grid.text(round(max_color,2), x=1, y=0.24, hjust = 1, vjust = 1, gp = gpar(cex = 0.8))

  popViewport()
  popViewport()
}

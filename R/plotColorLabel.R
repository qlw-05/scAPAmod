#' plotColorLabel function
#'
#' This function allows you to display gene name or PA information.
#'
#' @param color_labels display gene name or PA information
#'
#' @keywords plotColorLabel
#' @export
#'
plotColorLabel <-
function(color_labels){

  # pushViewport(viewport(yscale = c(1, length(color_labels)+1)))
  pushViewport(viewport(layout = grid.layout(length(color_labels),1)))
  for (i in 1:length(color_labels)) {
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    grid.text(color_labels[i], x=0, y=0.5, hjust = 0.5, vjust = 0, gp = gpar(cex = 0.8))

  # grid.rect(x=0, y=1, width = unit(1, "npc"), height = 1,
  #           default.units = "native", vjust = 0, hjust = 0,gp = gpar(fill = color_labels, col="black"))
    popViewport()
  }

  popViewport()
}

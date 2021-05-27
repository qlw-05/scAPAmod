#' plotComHeatmap function
#'
#' This function allows you to draw the PAs expression heat map of each component.
#'
#' @param mat a matrix
#' @param breaks set spacing
#' @param pal set color
#' @param fixed_max set maxinum value
#'
#' @keywords plotComHeatmap
#' @export
#'
plotComHeatmap <-
function(mat, breaks, pal, fixed_max){
  nr <- nrow(mat)
  pushViewport(viewport(layout = grid.layout(length(levels(factor(mat$component))),1)))
  for(i in seq_along(levels(factor(mat$component)))){
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    dat <- as.matrix(mat[which(mat$component == as.numeric(levels(factor(mat$component))[i])),-ncol(mat)])
    grid.raster(matrix(pal[as.numeric(cut2(dat, breaks, fixed_max = fixed_max))], nrow = nrow(dat)),
              interpolate = F, width =1, height = 1)
    popViewport()
  }
  popViewport()
  }

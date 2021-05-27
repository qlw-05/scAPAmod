#' plotHeatmap function
#'
#' This function allows you to draw PAs expression heat map.
#'
#' @param mat a matrix
#' @param breaks set spacing
#' @param pal set color
#' @param fixed_max set maximum
#'
#' @keywords plotHeatmap
#' @export
#'
plotHeatmap <-
function(mat, breaks, pal, fixed_max){
  nr <- nrow(mat)
  grid.raster(matrix(pal[as.numeric(cut2(mat, breaks, fixed_max = fixed_max))], nrow = nr),
              interpolate = F, width =1, height = 1)}

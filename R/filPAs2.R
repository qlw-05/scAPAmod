#' filPAs2 function
#'
#' This function allows you to filter PAs and cells.
#'
#' @param PAcounts the expression of PA (>=2)
#' @param gene gene corresponding to PAcounts
#' @param ftr the location of PA corresponding to PAcounts
#' @param gn proportion of filtered cells
#' @param cn proportion of filtered genes or PAs
#'
#' @keywords filPAs2
#' @export
#'
filPAs2 <-
function(PAcounts, gene, ftr, gn = 4, cn = 10){
  zeronum1 <- apply(PAcounts, 1, function(x){
    length(which(is.na(x)))
  })
  # ggplot() +
  #   geom_histogram(aes(x=zeronum1),bins = 50) +
  #   geom_vline(xintercept = ncol(non3UTR.data)/gn,color = "red")

  tofit2 <- which(zeronum1 < (ncol(PAcounts)/gn))
  PAcounts.tmp <- PAcounts[tofit2,,drop=FALSE]
  zeronum2 <- apply(PAcounts.tmp, 2, function(x){
    length(which(is.na(x)))
  })

  tofit1 <- which(zeronum2 <= (nrow(PAcounts.tmp)/cn))
  PAcounts.tmp <- PAcounts.tmp[,tofit1,drop=FALSE]
  #########################################################
  # filter annotation in the same time
  tUTR.gene <- gene[tofit2]
  PAftr <- ftr[tofit2]
  anno.tmp <- data.frame(gene = tUTR.gene, ftr = PAftr)
  return(list(fildata = PAcounts.tmp, anno = anno.tmp))
}

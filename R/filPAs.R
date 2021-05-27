#' filPAs function
#'
#' This function allows you to filter PAs and cells.
#'
#' @param PAcounts the expression of PA (>=2)
#' @param scAPAtrapds.tmp a PACdataset
#' @param gn proportion of filtered cells
#' @param cn proportion of filtered genes or PAs
#'
#' @keywords filPAs
#' @export
#'
filPAs <-
function(PAcounts, scAPAtrapds.tmp, gn = 4, cn = 10){
  # table(colSums(as.matrix(pair.data)[,,drop=FALSE]))
  zeronum1 <- apply(PAcounts, 1, function(x){
    length(which(x==0))
  })
  # ggplot() +
  #   geom_histogram(aes(x=zeronum1),bins = 50) +
  #   geom_vline(xintercept = ncol(non3UTR.data)/gn,color = "red")

  tofit2 <- which(zeronum1 < (ncol(PAcounts)/gn))
  PAcounts.tmp <- PAcounts[tofit2,,drop=FALSE]
  zeronum2 <- apply(PAcounts.tmp, 2, function(x){
    length(which(x==0))
  })

  tofit1 <- which(zeronum2 <= (nrow(PAcounts.tmp)/cn))
  PAcounts.tmp <- PAcounts.tmp[,tofit1,drop=FALSE]
  #########################################################
  # filter annotation in the same time
  tUTR.gene <- scAPAtrapds.tmp@anno$gene[tofit2]
  PAftr <- scAPAtrapds.tmp@anno$ftr[tofit2]
  PAchr <- scAPAtrapds.tmp@anno$chr[tofit2]
  PAstrand <- scAPAtrapds.tmp@anno$strand[tofit2]
  PAcoord <- scAPAtrapds.tmp@anno$coord[tofit2]
  PAftrStart <- scAPAtrapds.tmp@anno$ftr_start[tofit2]
  PAftrEnd <- scAPAtrapds.tmp@anno$ftr_end[tofit2]
  PAThreeL <- scAPAtrapds.tmp@anno$three_UTR_length[tofit2]
  anno.tmp <- data.frame(gene = tUTR.gene, ftr = PAftr, chr = PAchr, strand = PAstrand, coord = PAcoord,
                         ftr_start = PAftrStart, ftr_end = PAftrEnd, three_UTR_length = PAThreeL)
  return(list(fildata = PAcounts.tmp, anno = anno.tmp))
}

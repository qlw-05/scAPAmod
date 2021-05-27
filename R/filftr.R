#' filftr function
#'
#' This function allows you to filter out PAs in the intergenic region.
#'
#' @param scAPAtrapds a PACdataset
#'
#' @keywords filftr
#' @export
#'
filftr <-
function(scAPAtrapds){
  id <- which(scAPAtrapds@anno$ftr == "intergenic")
  pacanno <- cbind(as.character(scAPAtrapds@anno$chr[-id]) , as.character(scAPAtrapds@anno$strand[-id]) , scAPAtrapds@anno$coord[-id])
  colnames(pacanno) <- c("chr","strand","coord")
  pacfile <- cbind(pacanno,scAPAtrapds@counts[-id,])
  scPACpair <- readPACds(pacfile)
  scPACpair@colData <- scAPAtrapds@colData
  scPACpair@anno$gene <- scAPAtrapds@anno$gene[-id]
  scPACpair@anno$ftr <- scAPAtrapds@anno$ftr[-id]
  scPACpair@anno$ftr_start <- scAPAtrapds@anno$ftr_start[-id]
  scPACpair@anno$ftr_end <- scAPAtrapds@anno$ftr_end[-id]
  scPACpair@anno$three_UTR_length <- scAPAtrapds@anno$three_UTR_length[-id]
  return(scPACpair)
}

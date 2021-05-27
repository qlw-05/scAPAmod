#' filgene function
#'
#' This function allows you to remove genes named character(0).
#'
#' @param scAPAtrapds a PACdataset
#'
#' @keywords filgene
#' @export
#'
filgene <-
function(scAPAtrapds){
  id <- which(scAPAtrapds@anno$gene == "character(0)")
  pacanno <- cbind(as.character(scAPAtrapds@anno$chr[-id]) , as.character(scAPAtrapds@anno$strand[-id]) , scAPAtrapds@anno$coord[-id])
  colnames(pacanno) <- c("chr","strand","coord")
  pacfile <- cbind(pacanno,scAPAtrapds@counts[-id,])
  # pacfile$chr <- as.character(pacfile$chr)
  # pacfile$strand <- as.character(pacfile$strand)
  scPACpair <- readPACds(pacfile)
  scPACpair@colData <- data.frame(cell = scAPAtrapds@colData$Barcode, celltype = scAPAtrapds@colData$CellType)
  scPACpair@anno$gene <- scAPAtrapds@anno$gene[-id]
  scPACpair@anno$ftr <- scAPAtrapds@anno$ftr[-id]
  scPACpair@anno$ftr_start <- scAPAtrapds@anno$ftr_start[-id]
  scPACpair@anno$ftr_end <- scAPAtrapds@anno$ftr_end[-id]
  scPACpair@anno$three_UTR_length <- scAPAtrapds@anno$three_UTR_length[-id]
  return(scPACpair)
}

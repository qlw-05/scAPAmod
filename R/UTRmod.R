#' UTRmod function
#'
#' This function allows you to identifying modalities of APA usage in 3UTR.
#' The model used in the function is Gaussion Mixture model.
#'
#' @param APAdata a PACdataset
#' @param celltype the index of one same celltype
#'
#' @keywords UTRmod
#' @export
#'
UTRmod <-
function(APAdata,celltype){
  # extract pair PAs
  ex.data <- APAdata@counts
  all_anno <- APAdata@anno
  result <- extrPairPA(ex.data[,celltype], all_anno$gene, all_anno)
  # construction GMM model
  library(ClusterR)
  # msperm <- list(SCPUI,RSPUI,ESPUI)
  modalities <- getMod(PUIData = result$PUI)
  return(modalities)
}

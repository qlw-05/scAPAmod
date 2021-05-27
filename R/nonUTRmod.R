#' nonUTRmod function
#'
#' This function allows you to identifying modalities of APA usage in non 3UTR.
#' The model used in the function is Gaussion Mixture model.
#'
#' @param APAdata a PACdataset
#' @param celltype the index of one same celltype
#' @param gn proportion of filtered cells
#' @param cn proportion of filtered genes
#'
#' @keywords nonUTRmod
#' @export
#'
nonUTRmod <-
function(APAdata,celltype,gn = 4, cn = 10){
  # extract pair PAs
  ex.data <- APAdata@counts
  all_anno <- APAdata@anno
  result <- exnon3UTRPA(ex.data[,celltype], all_anno$gene,
                        all_anno, all_anno$ftr, gn, cn)
  # construction GMM model
  library(ClusterR)
  # msperm <- list(SCPUI,RSPUI,ESPUI)
  modalities <- getMod(PUIData = result$PUI)
  return(modalities)
}

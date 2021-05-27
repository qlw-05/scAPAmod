#' KStest function
#'
#' This function allows you to test Bimodal.
#'
#' @param PUIdata a PUI dataframe
#' @param tUTR.gene all gene in dataframe
#' @param ftr all ftr in dataframe
#' @param bigene gene which belong to Bimodal
#' @param bilabel label which belong to Bimodal
#'
#' @keywords KStest
#' @export
#'
KStest <- function(PUIdata, tUTR.gene, ftr, bigene, bilabel){

  p.val2 <- unlist(lapply(c(1:length(bigene)), function(x){
    if(!is.na(bigene[x])){
      pui <- as.matrix(PUIdata)[which(tUTR.gene %in% bigene[x]),]
      pui <- na.omit(pui[which(ftr[which(tUTR.gene %in% bigene[x])] != "3UTR"),])
      # if(length(which(colSums(pui) == 0)) > 0){
      #   pui <- as.matrix(pui[,-which(colSums(pui) == 0)])
      # }
      pval <- suppressWarnings(ks.test(pui[names(which(bilabel[[x]] == 0))], pui[names(which(bilabel[[x]] == 1))], exact = F))
      return(pval$p.value)
      }
    }))
  return(p.val2)
}

#' chisqtest function
#'
#' This function allows you to test Bimodal.
#'
#' @param filter.data a PA dataframe
#' @param tUTR.gene all gene in dataframe
#' @param bigene gene which belong to Bimodal
#' @param bilabel label which belong to Bimodal
#'
#' @keywords chisqtest
#' @export
#'
chisqtest <- function(filter.data, tUTR.gene, bigene, bilabel){
  p.val1 <- unlist(lapply(c(1:length(bigene)), function(x){

    pair <- as.matrix(na.omit(filter.data[which(tUTR.gene %in% bigene[x]),]))

    if(length(which(colSums(pair) == 0)) != ncol(pair)){
      l1 <- length(which(bilabel[[x]] == 0))
      l2 <- length(which(bilabel[[x]] == 1))
      if (which.max(c(l1,l2)) == 1) {
        for (y in 1:nrow(pair)) {
          # y = 1
          pair[y,1] <- sum(pair[y,names(which(bilabel[[x]] == 0))])*l2/l1
          pair[y,2] <- sum(pair[y,names(which(bilabel[[x]] == 1))])
        }
      }else{
        for (y in 1:nrow(pair)) {
          # y = 1
          pair[y,1] <- sum(pair[y,names(which(bilabel[[x]] == 0))])
          pair[y,2] <- sum(pair[y,names(which(bilabel[[x]] == 1))])*l1/l2
        }
      }
      pair <- pair[,1:2]
      colnames(pair) <- c(0,1)
      # print(pair)
      pval <- chisq.test(as.table(pair))
      return(pval$p.value)}}))

  return(p.val1)
}


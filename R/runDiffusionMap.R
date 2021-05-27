#' runDiffusionMap function
#'
#' This function allows you to sort by expression
#'
#' @param mat a Matrix
#'
#' @keywords runDiffusionMap
#' @export
#'
runDiffusionMap <-
function(mat){

  e <- try({
    dif <- DiffusionMap(as.ExpressionSet(as.data.frame(mat)));
    cat(sprintf("Eigenvalue of DC1: %f\n", eigenvalues(dif)[1]));
    res <- dif@eigenvectors[, "DC1"];
  }, silent=TRUE)

  if(class(e) == "try-error"){
    cat("There was a problem when running diffusion map. Trying PCA instead...\n")
    e <- try({
      pca =  prcomp(log10(mat+1),
                    center = TRUE,
                    scale. = TRUE);
      cat(sprintf("The standard deviations of PC1: %f\n", pca$sdev[1]));
      res <- pca$x[, "PC1"];
    }, silent=TRUE)

    if(class(e) == "try-error"){
      res <- rep(1, nrow(mat))
    }else{
      return(res)
    }
  }else{
    return(res)
  }
}

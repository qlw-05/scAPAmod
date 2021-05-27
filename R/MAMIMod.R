#' MAMIMod function
#'
#' This function allows you to analyze APA's preferences.
#' The model used in the function is Gaussion Mixture model.
#'
#' @param APAdata a PACdataset
#' @param celltype the index of one same celltype
#' @param type select "MajorPA" or "MinorPA"
#'
#' @keywords MAMIMod
#' @export
#'
MAMIMod <-
function(APAdata,celltype,type){
  result <- exMajorPA(APAdata, celltype)
  if(type == "MajorPA"){
    PAmax <- result$PAmax
    modalities <- getMMod(PAmax,"PAmax")}
  if(type == "MinorPA"){
    PAmin <- result$PAmin
    modalities <- getMMod(PAmin,"PAmin")
  }
  return(modalities)
}

#' exMajorPA function
#'
#' This function allows you to calculate the ratio of PAs and select the maximum and minimum values.
#'
#' @param scAPAtrapds a PACdataset
#' @param cell_type select the cell type
#'
#' @keywords exMajorPA
#' @export
#'
exMajorPA <-
function(scAPAtrapds, cell_type = "SC"){
  scAPAtrapds.tmp <- scAPAtrapds
  if (length(which(scAPAtrapds.tmp@anno$gene == "character(0)")) > 0 ) {
    scAPAtrapds.tmp <- filgene(scAPAtrapds.tmp)
  }
  if(length(which(scAPAtrapds.tmp@anno$ftr == "intergenic")) > 0 ){
    scAPAtrapds.tmp <- filftr(scAPAtrapds.tmp)
  }
  PAcounts <- scAPAtrapds.tmp@counts[,which(scAPAtrapds.tmp@colData$celltype == cell_type)]

  # extract max and min PAs
  unigene <- unique(scAPAtrapds.tmp@anno$gene)
  id <- match(scAPAtrapds.tmp@anno$gene, unigene)
  unigene2 <- unigene[which(table(id) >= 3)]
  ## don't divide 3UTR and non3UTR
  mmPAs <- lapply(unigene2, function(g){
    # g<- unigene2[404]
    index <- which(scAPAtrapds.tmp@anno$gene %in% g)
    PAs <- PAcounts[index,]
    ftrs <- scAPAtrapds.tmp@anno$ftr[index]
    PAratio <- apply(PAs, 2, function(p){p/sum(p)})
    maxID <- which.max(rowSums(PAs))
    minID <- which.min(rowSums(PAs))

    maxPA <- as.matrix(PAratio[maxID,])
    minPA <- as.matrix(PAratio[minID,])
    majorftr <- as.character(ftrs[maxID])
    minorftr <- as.character(ftrs[minID])

    return(list(maxPA = maxPA, minPA = minPA, majorftr = majorftr, minorftr = minorftr,gene = g))
  })

  maxPAs <- lapply(mmPAs, function(m){
    return(t(m$maxPA))
  })
  maxPAs <- do.call("rbind", maxPAs)
  maxftrs <- unlist(lapply(mmPAs, function(m){return(m$majorftr)}))
  maxgene <- unlist(lapply(mmPAs, function(m){m$gene}))
  PAmax <- filPAs2(maxPAs, maxgene, ftr = maxftrs, gn = 4)

  minPAs <- lapply(mmPAs, function(m){
    return(t(m$minPA))
  })
  minPAs <- do.call("rbind", minPAs)
  minftrs <- unlist(lapply(mmPAs, function(m){return(m$minorftr)}))
  mingene <- unlist(lapply(mmPAs, function(m){return(m$gene)}))
  PAmin <- filPAs2(minPAs, mingene, ftr = minftrs, gn = 4)
  return(list(PAmax = PAmax, PAmin = PAmin))
}

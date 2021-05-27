#' extrPairPA function
#'
#' This function allows you to filter PAs and cells, and calculate the PUI of genes.
#'
#' @param exc1.g1 counts from a PACdataset(APA expression in the same cell type)
#' @param all_gene gene from anno in PACdataset
#' @param all_anno anno from a PACdataset
#' @param gn proportion of filtered cells
#' @param cn proportion of filtered genes
#'
#' @keywords extrPairPA
#' @export
#'
extrPairPA <-
function(exc1.g1, all_gene, all_anno, gn=4, cn=10){
  uni.gene <- sort(unique(all_gene))
  index <- match(all_gene,uni.gene)
  ########################################################
  # pair data
  snum <- table(index)
  pair.gene <- uni.gene[which(snum>=2)]
  pair.data <- exc1.g1[which(all_gene %in% pair.gene),]
  pair.data.gene <- all_gene[which(all_gene %in% pair.gene)]
  #########################################################
  # add annotation information(constract PACdataSet object)
  PAchr <- all_anno$chr[which(all_gene %in% pair.gene)]
  PAstrand <- all_anno$strand[which(all_gene %in% pair.gene)]
  PAcoord <- all_anno$coord[which(all_gene %in% pair.gene)]
  PAftrStart <- all_anno$ftr_start[which(all_gene %in% pair.gene)]
  PAftrEnd <- all_anno$ftr_end[which(all_gene %in% pair.gene)]
  PAThreeL <- all_anno$three_UTR_length[which(all_gene %in% pair.gene)]
  #########################################################
  # filter low expression cells
  # table(colSums(as.matrix(pair.data)[,,drop=FALSE]))
  zeronum1 <- apply(pair.data, 1, function(x){
    length(which(x==0))
  })

  tofit2 <- which(zeronum1 < (ncol(pair.data)/gn))
  tUTR.pair.cd <- pair.data[tofit2,,drop=FALSE]
  zeronum2 <- apply(tUTR.pair.cd, 2, function(x){
    length(which(x==0))
  })

  tofit1 <- which(zeronum2 <= (nrow(tUTR.pair.cd)/cn))
  tUTR.pair.cd.tmp <- tUTR.pair.cd[,tofit1,drop=FALSE]
  #########################################################
  # filter annotation in the same time
  tUTR.gene <- pair.data.gene[tofit2]
  PAchr <- PAchr[tofit2]
  PAstrand <- PAstrand[tofit2]
  PAcoord <- PAcoord[tofit2]
  PAftrStart <- PAftrStart[tofit2]
  PAftrEnd <- PAftrEnd[tofit2]
  PAThreeL <- PAThreeL[tofit2]
  ########################################################
  # trimmed single PA
  geneid <- unique(tUTR.gene)
  index <- match(tUTR.gene,geneid)
  snum1 <- table(index)
  pair.gene1 <- geneid[which(snum1>=2)]

  tUTR.pair.cd.tmp <- tUTR.pair.cd.tmp[which(tUTR.gene %in% pair.gene1),]
  tUTR.gene <- tUTR.gene[which(tUTR.gene %in% pair.gene1)]

  PAchr <- PAchr[which(tUTR.gene %in% pair.gene1)]
  PAstrand <- PAstrand[which(tUTR.gene %in% pair.gene1)]
  PAcoord <- PAcoord[which(tUTR.gene %in% pair.gene1)]
  PAftrStart <- PAftrStart[which(tUTR.gene %in% pair.gene1)]
  PAftrEnd <- PAftrEnd[which(tUTR.gene %in% pair.gene1)]
  PAThreeL <- PAThreeL[which(tUTR.gene %in% pair.gene1)]
  ########################################################
  # trimming character(0) of gene
  if(length(which(tUTR.gene == "character(0)")) > 0){
    tUTR.pair.cd.tmp <- tUTR.pair.cd.tmp[-which(tUTR.gene == "character(0)"),]
    PAchr <- PAchr[-which(tUTR.gene == "character(0)")]
    PAstrand <- PAstrand[-which(tUTR.gene == "character(0)")]
    PAcoord <- PAcoord[-which(tUTR.gene == "character(0)")]
    PAUPAstart <- PAUPAstart[-which(tUTR.gene == "character(0)")]
    PAUPAend <- PAUPAend[-which(tUTR.gene == "character(0)")]
    PAftrStart <- PAftrStart[-which(tUTR.gene == "character(0)")]
    PAftrEnd <- PAftrEnd[-which(tUTR.gene == "character(0)")]
    PAThreeL <- PAThreeL[-which(tUTR.gene == "character(0)")]
    tUTR.gene <- tUTR.gene[-which(tUTR.gene == "character(0)")]
  }
  ########################################################
  ##calculate PUI
  pacanno <- cbind(as.character(PAchr), as.character(PAstrand), as.integer(PAcoord))
  colnames(pacanno) <- c("chr","strand","coord")
  pacfile <- cbind(pacanno,tUTR.pair.cd.tmp)
  pacfile$chr <- as.character(pacfile$chr)
  pacfile$strand <- as.character(pacfile$strand)
  scPACpair <- readPACds(pacfile)
  scPACpair@anno$gene <- tUTR.gene
  scPACpair@anno$coord <- as.integer(as.character(scPACpair@anno[["coord"]]))
  scPACpair@anno$ftr <- as.character(rep("3UTR",nrow(tUTR.pair.cd.tmp)))
  scPACpair@anno$ftr_start <- as.integer(PAftrStart)
  scPACpair@anno$ftr_end <- as.integer(PAftrEnd)
  scPACpair@anno$three_UTR_length <- PAThreeL

  # calculate every PA PUI
  # PUI <- movPAindex(scPACpair, method = "geo")
  # calculate every gene PUI
  PUI <- movAPAindex(scPACpair, method = "GPI")

  return(list(PUI = PUI, filter.data = tUTR.pair.cd.tmp, gene = tUTR.gene))
}

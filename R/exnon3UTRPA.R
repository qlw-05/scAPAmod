#' exnon3UTRPA function
#'
#' This function allows you to filter PAs and cells, and calculate the PUI of PAs.
#'
#' @param exc1.g1 counts from a PACdataset(APA expression in the same cell type)
#' @param all_gene gene from anno in PACdataset
#' @param all_anno anno from a PACdataset
#' @param all_ftr ftr from anno in PACdataset
#' @param gn proportion of filtered cells
#' @param cn proportion of filtered genes
#'
#' @keywords exnon3UTRPA
#' @export
#'
exnon3UTRPA <-
function(exc1.g1, all_gene, all_anno, all_ftr, gn=4, cn=10){
  uni.gene <- sort(unique(all_gene))
  index <- match(all_gene, uni.gene)
  ########################################################
  # pair data
  snum <- table(index)
  pair.gene <- uni.gene[which(snum>=2)]
  #########################################################
  #grab non3UTR which includes 3UTR
  non3UTR.data <- data.frame()
  non3UTR.gene <- c()
  non3UTR.ftr <- c()
  PAchrs <- c()
  PAstrands <- c()
  PAcoords <- c()
  PAftrStarts <- c()
  PAftrEnds <- c()
  PAThreeLs <- c()
  non3UTRinfo <- lapply(pair.gene, function(m){
    # m <- pair.gene[4]
    pair.data.ftr <- all_ftr[which(all_gene %in% m)]
    if(length(which(pair.data.ftr == "3UTR")) > 0){
      pair.data <- exc1.g1[which(all_gene %in% m),]
      pair.data1 <- pair.data[-which(pair.data.ftr == "3UTR"),]
      #=========================================================#
      chr <- all_anno$chr[which(all_gene %in% m)]
      strand <- all_anno$strand[which(all_gene %in% m)]
      coord <- all_anno$coord[which(all_gene %in% m)]
      ftrStart <- all_anno$ftr_start[which(all_gene %in% m)]
      ftrEnd <- all_anno$ftr_end[which(all_gene %in% m)]
      ThreeL <- all_anno$three_UTR_length[which(all_gene %in% m)]
      #=========================================================#

      if(nrow(pair.data1) != 0){
        #=========================================================#
        PAchr <- chr[-which(pair.data.ftr == "3UTR")]
        PAstrand <- strand[-which(pair.data.ftr == "3UTR")]
        PAcoord <- coord[-which(pair.data.ftr == "3UTR")]
        PAftrStart <- ftrStart[-which(pair.data.ftr == "3UTR")]
        PAftrEnd <- ftrEnd[-which(pair.data.ftr == "3UTR")]
        PAThreeL <- ThreeL[-which(pair.data.ftr == "3UTR")]
        #=========================================================#
        non3UTR <- pair.data1[which.max(rowSums(pair.data1)),]
        non3UTR <- rbind(pair.data[which(pair.data.ftr == "3UTR"),],non3UTR)
        non3UTR.gene <- c(non3UTR.gene, rep(m, nrow(non3UTR)))
        non3UTR.data <- rbind(non3UTR.data, non3UTR)
        ftr <- pair.data.ftr[-which(pair.data.ftr == "3UTR")][which.max(rowSums(pair.data1))]
        non3UTR.ftr <- c(non3UTR.ftr, c(rep("3UTR", length(which(pair.data.ftr == "3UTR"))), ftr))
        PAchrs <- c(PAchrs, c(as.character(chr[which(pair.data.ftr == "3UTR")]), as.character(PAchr[which.max(rowSums(pair.data1))])))
        PAstrands <- c(PAstrands, c(as.character(strand[which(pair.data.ftr == "3UTR")]), as.character(PAstrand[which.max(rowSums(pair.data1))])))
        PAcoords <- c(PAcoords, c(coord[which(pair.data.ftr == "3UTR")], PAcoord[which.max(rowSums(pair.data1))]))
        PAftrStarts <- c(PAftrStarts, c(ftrStart[which(pair.data.ftr == "3UTR")], PAftrStart[which.max(rowSums(pair.data1))]))
        PAftrEnds <- c(PAftrEnds, c(ftrEnd[which(pair.data.ftr == "3UTR")], PAftrEnd[which.max(rowSums(pair.data1))]))
        PAThreeLs <- c(PAThreeLs, c(ThreeL[which(pair.data.ftr == "3UTR")], PAThreeL[which.max(rowSums(pair.data1))]))

      }
    }
    return(list(non3UTR.data = non3UTR.data, non3UTR.gene = non3UTR.gene, non3UTR.ftr = non3UTR.ftr, PAchr = PAchrs,
                PAstrand = PAstrands, PAcoord = PAcoords, PAftrStart = PAftrStarts, PAftrEnd = PAftrEnds, PAThreeL = PAThreeLs))
  })
  non3UTR.data <- do.call("rbind", lapply(non3UTRinfo, function(x){x$non3UTR.data}))
  pair.data.gene <- unlist(lapply(non3UTRinfo, function(y){y$non3UTR.gene}))
  PAftr <- unlist(lapply(non3UTRinfo, function(y){y$non3UTR.ftr}))
  PAchr <- unlist(lapply(non3UTRinfo, function(y){y$PAchr}))
  PAstrand <- unlist(lapply(non3UTRinfo, function(y){y$PAstrand}))
  PAcoord <- unlist(lapply(non3UTRinfo, function(y){y$PAcoord}))
  PAftrStart <- unlist(lapply(non3UTRinfo, function(y){y$PAftrStart}))
  PAftrEnd <- unlist(lapply(non3UTRinfo, function(y){y$PAftrEnd}))
  PAThreeL <- unlist(lapply(non3UTRinfo, function(y){y$PAThreeL}))
  #########################################################
  # filter low expression cells
  # table(colSums(as.matrix(pair.data)[,,drop=FALSE]))
  zeronum1 <- apply(non3UTR.data, 1, function(x){
    length(which(x==0))
  })
  # ggplot() +
  #   geom_histogram(aes(x=zeronum1),bins = 50) +
  #   geom_vline(xintercept = ncol(non3UTR.data)/gn,color = "red")

  tofit2 <- which(zeronum1 < (ncol(non3UTR.data)/gn))
  non3UTR.data.tmp <- non3UTR.data[tofit2,,drop=FALSE]
  zeronum2 <- apply(non3UTR.data.tmp, 2, function(x){
    length(which(x==0))
  })

  tofit1 <- which(zeronum2 <= (nrow(non3UTR.data.tmp)/cn))
  non3UTR.data.tmp <- non3UTR.data.tmp[,tofit1,drop=FALSE]
  #########################################################
  # filter annotation in the same time
  tUTR.gene <- pair.data.gene[tofit2]
  PAftr <- PAftr[tofit2]
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
  if(length(which(pair.gene1 == "character(0)")) > 1){
      pair.gene1 <- pair.gene1[-which(pair.gene1 == "character(0)")]
  }
  filnon3UTR <- lapply(pair.gene1, function(g){
    # g <- pair.gene1[2]

    if(length(unique(PAftr[which(tUTR.gene %in% g)])) == 2){
      non3UTR.fildata <- non3UTR.data.tmp[which(tUTR.gene %in% g),]
      fil.gene <- tUTR.gene[which(tUTR.gene %in% g)]
      PA.ftr <- PAftr[which(tUTR.gene %in% g)]
      PA.chr <- PAchr[which(tUTR.gene %in% g)]
      PA.strand <- PAstrand[which(tUTR.gene %in% g)]
      PA.coord <- PAcoord[which(tUTR.gene %in% g)]
      PA.ftrStart <- PAftrStart[which(tUTR.gene %in% g)]
      PA.ftrEnd <- PAftrEnd[which(tUTR.gene %in% g)]
      PA.ThreeL <- PAThreeL[which(tUTR.gene %in% g)]
      return(list(non3UTR.data = non3UTR.fildata, non3UTR.gene = fil.gene, non3UTR.ftr = PA.ftr, PAchr = PA.chr,
                  PAstrand = PA.strand, PAcoord = PA.coord, PAftrStart = PA.ftrStart, PAftrEnd = PA.ftrEnd, PAThreeL = PA.ThreeL))
    }

  })
  non3UTR.fildata <- do.call("rbind", lapply(filnon3UTR, function(x){x$non3UTR.data}))
  tUTR.gene.fil <- unlist(lapply(filnon3UTR, function(y){y$non3UTR.gene}))
  PAftrs <- unlist(lapply(filnon3UTR, function(y){y$non3UTR.ftr}))
  PAchrs <- unlist(lapply(filnon3UTR, function(y){y$PAchr}))
  PAstrands <- unlist(lapply(filnon3UTR, function(y){y$PAstrand}))
  PAcoords <- unlist(lapply(filnon3UTR, function(y){y$PAcoord}))
  PAftrStarts <- unlist(lapply(filnon3UTR, function(y){y$PAftrStart}))
  PAftrEnds <- unlist(lapply(filnon3UTR, function(y){y$PAftrEnd}))
  PAThreeLs <- unlist(lapply(filnon3UTR, function(y){y$PAThreeL}))

  ########################################################
  ##calculate PUI
  pacanno <- cbind(as.character(PAchrs), as.character(PAstrands), as.integer(PAcoords))
  colnames(pacanno) <- c("chr","strand","coord")
  pacfile <- cbind(pacanno,non3UTR.fildata)
  pacfile$chr <- as.character(pacfile$chr)
  pacfile$strand <- as.character(pacfile$strand)
  scPACpair <- readPACds(pacfile)
  scPACpair@anno$gene <- tUTR.gene.fil
  scPACpair@anno$coord <- as.integer(as.character(scPACpair@anno[["coord"]]))
  scPACpair@anno$ftr <- as.character(PAftrs)
  scPACpair@anno$ftr_start <- as.integer(PAftrStarts)
  scPACpair@anno$ftr_end <- as.integer(PAftrEnds)
  scPACpair@anno$three_UTR_length <- PAThreeLs

  # calculate every PA PUI
  PUI <- movPAindex(scPACpair, method = "geo")
  # calculate every gene PUI
  # PUI <- movAPAindex(scPACpair, method = "GPI")
  if (nrow(PUI) != length(tUTR.gene.fil)) {
    t <- tUTR.gene.fil[which(PAftrs == "intergenic")]
    PUI.gene <- tUTR.gene.fil[-which(tUTR.gene.fil %in% t)]
    ftrs <- PAftrs[-which(tUTR.gene.fil %in% t)]
    non3UTR.fildata <- non3UTR.fildata[-which(tUTR.gene.fil %in% t),]
  }

  return(list(PUI = PUI, filter.data = non3UTR.fildata, gene = PUI.gene, ftr = ftrs))
}

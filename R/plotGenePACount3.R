#' plotGenePACount3 function
#'
#' This function allows you to draw PAs expression heat map and expression heat map of each component.
#'
#' @param org.Mm.eg.db select an annotation information for a species
#' @param id index belonging to Bimodal(ENTREZID)
#' @param tUTR.pair.cd.tmp PAs expression(filter.data-output of extrPairPA)
#' @param tUTR.gene.fil gene sample
#' @param label label of cell components
#'
#' @keywords plotGenePACount3
#' @export
#'
library(grid)
plotGenePACount3 <-
function(org.Mm.eg.db, id, tUTR.pair.cd.tmp, tUTR.gene.fil,label){
  heights <- c(12,2)
  col_ratio <- unit(c(1,1,4,1), c("null","cm", "null", "null"))
  breaks <- 100
  pal <- colorRampPalette(c('white','blue'))
  pal <- pal(breaks)

  same.gene <- select(org.Mm.eg.db, keys = id, columns = c("SYMBOL","ENSEMBL","GENENAME"),
                      keytype = "ENTREZID")
  ex.samegene <- as.data.frame(t(tUTR.pair.cd.tmp[which(tUTR.gene.fil %in% id),]))
  if (length(which(rowSums(ex.samegene) == 0) ) > 0) {
    ex.samegene <- ex.samegene[-which(rowSums(ex.samegene) == 0),]
  }
  if(length(is.na(ex.samegene[,2])) > 0 ){ex.samegene <- na.omit(ex.samegene)}
  ex.samegene$component <- label[match(rownames(ex.samegene), names(label))]
  # if(length(which(rowSums(ex.samegene) == 0)) > 0){
  #   ex.samegene <- ex.samegene[-which(colSums(ex.samegene) == 0),]
  # }
  sel_danger <- rowSums(ex.samegene[,-ncol(ex.samegene)]) == 0
  sel_safe <- !(sel_danger)
  sel_safe[which(sel_danger)[1]] <- TRUE
  dc1_value <- runDiffusionMap(ex.samegene[sel_safe,-ncol(ex.samegene)])
  dc1 <- numeric(nrow(ex.samegene))
  dc1[sel_safe] <- dc1_value
  dc1[!sel_safe] <- dc1[which(sel_danger)[1]]
  neworder <- order(dc1, decreasing = T)
  ex.samegene <- ex.samegene[neworder,]
  # colorlabel <- colorlabel[neworder,]
  ex.all <- as.matrix(ex.samegene[,-ncol(ex.samegene)])


  grid.newpage()
  pushViewport(plotViewport(c(1,0,1,0)))
  # design plot
  pushViewport(viewport(layout = grid.layout(3, 1, heights = heights)))

  # first plot
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  pushViewport(plotViewport(c(0.5,0,0.5,0)))
  pushViewport(viewport(layout = grid.layout(2, length(col_ratio),
                                             heights = unit(c(2,12), c("lines", "null") ),
                                             widths=col_ratio)))

  # plot PA heatmap
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
  plotHeatmap(ex.all ,breaks, pal, max(ex.all))
  popViewport()


  # add a label about expressions of PAs
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
  plotColorKey2(pal[1:breaks], min(ex.all), max(ex.all))
  popViewport()

  # add PA info label
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  plotColorLabel(rev(unique(same.gene$SYMBOL)))
  # plotColorLabel(rev(colorlabel))
  popViewport()

  popViewport()
  popViewport()
  popViewport()

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  # pushViewport(viewport(layout = grid.layout(1, length(col_ratio), heights = unit(2,"lines"), widths=col_ratio)))
  # pushViewport(viewport(layout = grid.layout(2, length(col_ratio), widths=col_ratio)))
  pushViewport(viewport(layout = grid.layout(2, length(col_ratio),
                                             heights = unit(c(2,5), "lines"),
                                             widths=col_ratio )))

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
  # plotHeatmap(as.matrix(ex.samegene), breaks, pal, max(ex.samegene[,-ncol(ex.samegene)]))
  plotComHeatmap(ex.samegene,breaks, pal, max(ex.samegene[,-ncol(ex.samegene)]))
  popViewport()


  # add a label about expressions of PAs
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 4))
  com <- levels(as.factor(label))
  x.min <- c(NULL)
  x.max <- c(NULL)
  for (m in 1:length(com)) {
    mt <- ex.all[which(ex.samegene[,ncol(ex.samegene)] == com[m]),]
    x.min <- c(x.min, min(mt))
    x.max <- c(x.max, max(mt))
  }
  # plotYvalue2(length(com),pal[1:breaks], x.max, x.min)
  plotYvalue(length(com), x.max, x.min)
  popViewport()

  # add PA info label
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  plotColorLabel(com)
  popViewport()

  popViewport()
  popViewport()


  popViewport()
  popViewport()

}

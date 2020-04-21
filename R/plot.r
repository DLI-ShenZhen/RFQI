#' plot EIC
#' 
#' @param f.in path of mzxml files
#' @param refFile path of reference mzXML files, using for retention time alignment
#' @param mzr m/z range, if appointed, plot mzr area EIC
#' @param rtr rt range, if appointed, plot rtr area EIC
#' @param features results from function getPeaks_multipleFiles, containing two data frame: raw peaks and features
#' @param idx id of features from featureTable, chrom_plot only plot idx corresponding EIC
#' @param col colours for each sample, if not appointed, using brewer.pal(9, "Set3")
#' @return a plot
#' @export

chrom_plot = function (f.in, refFile = NULL, features, idx, mzr = NULL, rtr = NULL, col=NULL){
  require(xcms)
  require(RColorBrewer)
  if (is.null(col)){colors = brewer.pal(9, "Set3")} else {colors=col}
  MSData = list()
  for (i in 1:length(f.in)) {
    if (!is.null(refFile)) {
      f.in.local = c(refFile, f.in[i])
      raw_data = readMSData(files = f.in.local, mode = "onDisk", 
                            msLevel. = NULL, centroided. = TRUE)
      rtcor = adjustRtime(raw_data, param = ObiwarpParam(binSize = 0.25, 
                                                         centerSample = 1, localAlignment = FALSE), msLevel = 1L)
      Biobase::featureData(raw_data)@data$retentionTime = rtcor
    }
    else {
      raw_data = readMSData(files = f.in[i], mode = "onDisk", 
                            msLevel. = NULL, centroided. = TRUE)
    }
    MSData[[length(MSData) + 1]] = raw_data
  }
  for (id in idx) {
    id_feature = features$group_FT[id, , drop = TRUE]
    id_peakId = unlist(strsplit(id_feature$id, split = ","))
    id_peakId = as.numeric(id_peakId)
    id_peaks = features$peaks[id_peakId, , drop = FALSE]
    if (is.null(mzr)) 
      mzr = as.numeric(id_feature[c("mzmin", "mzmax")])
    if (is.null(rtr)) 
      rtr = as.numeric(id_feature[c("rtmin", "rtmax")])
    if (diff(mzr) == 0) {
      next
    }
    dataMat = matrix(ncol = 3)
    colnames(dataMat) = c("rtime", "intensity", "color")
    fileIdx = ifelse(is.null(refFile), 1, 2)
    for (i in 1:length(f.in)) {
      chr_raw = chromatogram(MSData[[i]], mz = mzr, rt = rtr, 
                             aggregationFun = "sum")
      rtime = chr_raw@.Data[[fileIdx]]@rtime
      intensity = chr_raw@.Data[[fileIdx]]@intensity
      dataMat = rbind(dataMat, data.frame(rtime = rtime, 
                                          intensity = intensity, color = colors[i]))
    }
    dataMat = dataMat[-1, ]
    xlim = range(dataMat[, 1], na.rm = TRUE)
    ylim = range(dataMat[, 2], na.rm = TRUE)
    if (Inf %in% ylim) {
      next
    }
    ylim[1] = 0
    for (i in 1:length(f.in)) {
      if (i == 1) {
        plot(dataMat[dataMat[, "color"] == colors[i], 
                     1], dataMat[dataMat[, "color"] == colors[i], 
                                 2], type = "l", col = colors[i], xlim = xlim, 
             ylim = ylim, lwd = 2, xlab = "rtime", ylab = "intensity")
      }
      else {
        lines(dataMat[dataMat[, "color"] == colors[i], 
                      1], dataMat[dataMat[, "color"] == colors[i], 
                                  2], type = "l", col = colors[i], lwd = 2)
      }
    }
    # title(main = sprintf("mz: %s ~ %s", round(mzr[1], 3), 
    #                      round(mzr[2], 3)))
  }
  
  return (dataMat)
}

#' plot MS2 Cluster
#' 
#' plot MS2 similarity network using igraph
#' 
#' @param df correlation matrix
#' @param cutoff similarity score threshold
#' @param main main title of plot
#' @param lable if TRUE, add similarity score on edges
#' @export
plotMS2Cluster = function(df, cutoff=0.4, main="", label=TRUE){
  
  len = nrow(df) - 1
  NodeEdge = lapply(1:len, function(i){
    x = df[i,(i+1):nrow(df)]
    from = rep(i,length(x))
    to = (i+1):nrow(df)
    cbind(from=from, to=to, weight=x)
  })
  NodeEdge = do.call(rbind,NodeEdge)
  NodeEdge = NodeEdge[NodeEdge[,"weight"]>cutoff,]
  
  g = graph_from_data_frame(NodeEdge, directed = FALSE)
  
  if (label){
    plot.igraph(g, vertex.size=1, shape="circle",
                label.font=2, label.cex=15, label.dist=1, label.color="black", label.family="bold",
                edge.label=round(NodeEdge[,"weight"],2),
                edge.label.color = "red",
                edge.width = NodeEdge[,"weight"],
                edge.label.cex = 0.8,
                layout=layout_nicely,
                main=main)
  } else {
    plot.igraph(g, vertex.size=1, shape="circle",
                label.font=2, label.cex=15, label.dist=1, label.color="black", label.family="bold",
                layout=layout_nicely,
                main=main)
  }
  
}



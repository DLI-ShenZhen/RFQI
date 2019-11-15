#' summary_L is to compute rt and ppm shift
#' 
#' @param df data frame whose colnames containing c("mzmin","mzmax","rtmin","rtmax")
#' @return numeric matrix containing m/z difference, rt difference, and m/z difference in ppm unit
#' @export
summary_L = function(df){
  df = df[,c("mzmin","mzmax","rtmin","rtmax")]
  df = as.matrix(df); mode(df) = "numeric"
  mz.diff = df[,"mzmax"] - df[,"mzmin"]
  rt.diff = df[,"rtmax"] - df[,"rtmin"]
  mz.diff.ppm = mz.diff / df[,"mzmax"] *10e6
  return(cbind(mz.diff=mz.diff, rt.diff=rt.diff, mz.diff.ppm))
}

#' Function plotArea is used to plot feature area of one FT or two
#' 
#' @importFrom scales alpha
#' @param refPeaks dataframe or matrix containing columns of c("mzmin","mzmax","rtmin","rtmax")
#' @param group factor, length as row of refPeaks, indicating every features' group
#' @param main title of plot
#' @param mzr m/z range
#' @param rtr rt range
#' @return plot
#' @export
plotArea = function(refPeaks, group=NULL, main="", mzr=NULL, rtr=NULL){
  # refPeaks = features$group_FT
  
  area = refPeaks[,c("mzmin", "mzmax", "rtmin", "rtmax")]
  area = as.matrix(area)
  mode(area) = "numeric"
  
  xlim = range(unlist(area[,c("rtmin","rtmax")]))
  ylim = range(unlist(area[,c("mzmin","mzmax")]))
  
  # mzr and rtr filter
  if (! is.null(mzr) | !is.null(rtr)){
    Mz = ifelse(is.null(mzr), FALSE, TRUE)
    Rt = ifelse(is.null(rtr), FALSE, TRUE)
    
    # reflesh xlim and ylim
    if(is.null(mzr)){mzr=ylim}
    if(is.null(rtr)){rtr=xlim}
    ref = c(mzr, rtr)
    
    xlim = rtr
    ylim = mzr
    
    # filtr area
    flag = apply(area, 1, function(x){
      judgePeakOverlap(ref=ref, alt=x, Mz=Mz, Rt=Rt, absMz=0, absRt=0 )
    })
    area = area[flag,,drop=FALSE]
    
  }
  
  
  # plot(x=0,y=0,xlim=xlim, ylim=ylim, xlab="retention time", ylab="m/z", cex.axis=1, cex.lab=2, main="features region overlap", xaxt="n", yaxt="n")
  # axis(side=1, at <- seq(ceiling(xlim)[1], ceiling(xlim)[2],10), labels = at)
  # axis(side=2, at <- seq(ceiling(ylim)[1], ceiling(ylim)[2],5), labels = at)
  plot(x=0,y=0,xlim=xlim, ylim=ylim, xlab="retention time", ylab="m/z", cex.axis=1, cex.lab=1.5, main=main)
  
  # --- single group plot ----- #
  if(is.null(group)){
    # plot.window(xlim=xlim, ylim=ylim)
    rect(xleft=area[,"rtmin"], ybottom=area[,"mzmin"], xright=area[,"rtmax"], ytop=area[,"mzmax"], density=-1, angle=45, border=alpha("blue",0.5), lwd=0.5)
  }
  
  # --- two group plot -------- #
  if(!is.null(group)){
    group = as.factor(group)
    area.split = split.data.frame(area, f = group, drop=FALSE)
    a = area.split[[1]]
    b = area.split[[2]]
    rect(xleft=a[,"rtmin"], ybottom=a[,"mzmin"], xright=a[,"rtmax"], ytop=a[,"mzmax"], density=-1, angle=45, border=alpha("blue",0.5), lwd=1)
    rect(xleft=b[,"rtmin"], ybottom=b[,"mzmin"], xright=b[,"rtmax"], ytop=b[,"mzmax"], density=-1, angle=45, border=alpha("red",0.5), lwd=1)
    
  }
  
}

#' Function plotFeatureMap is to plot Feature, peaks, and MS2 point in a figure
#' 
#' @import pheatmap
#' @import grDevices 
#' @import graphics
#' @param idx index
#' @param MS2_inner_cor MS2 similarity matrix for MS2 clustering
#' @param features containing peaks and features
#' @param MS2DB db.MS2 object 
#' @param heatmap logical, if true, plot MS2 similarity heatmap, else plot featureMap
#' @return figure 
#' @export
plotFeatureMap = function(idx, MS2_inner_cor, features, MS2DB, heatmap=FALSE){
  db.MS2 = MS2DB
  
  idx_ms2Cor = MS2_inner_cor[[idx]]
  idx_MS2 = checkMS2(db.MS2, idx = idx)
  idx_Peaks = checkFT(features = features, idx=idx)
  
  cluster = get_ms2cluster(idx=idx, MS2_inner_cor = MS2_inner_cor, eps=0.5, minPts=2)
  k = length(cluster$cluster)
  col = colorRampPalette(c("blue","green"))(k)
  idx_MS2_clust = rep(c(1:k,0), c(sapply(cluster$cluster, length), length(cluster$noise)))
  names(idx_MS2_clust) = unlist(cluster)
  
  # plot feature map
  mzr = features$group_FT[idx,c("mzmin","mzmax"),drop=TRUE]; mode(mzr)="numeric"
  rtr = features$group_FT[idx,c("rtmin","rtmax"),drop=TRUE]; mode(rtr)="numeric"
  
  
  refPeaks = rbind(idx_Peaks[,c("mzmin","mzmax","rtmin","rtmax")], 
                   features$group_FT[idx,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE])
  group = c(rep("P",nrow(idx_Peaks)),"F")  
  
  plotArea(refPeaks = refPeaks,
           group = group,
           main = sprintf("%s\nmz.diff = %sDa (%sppm); rt.diff = %ss;", 
                          idx, 
                          round(diff(mzr), 3),
                          round(diff(mzr)/mzr[1]*10e6,0),
                          round(diff(rtr),0))
  )
  
  colors = col[idx_MS2_clust[as.character(1:length(idx_MS2_clust))]]
  points(idx_MS2$precursor[,"retentionTime"], idx_MS2$precursor[,"precursorMZ"], pch=16, cex=1, col=colors)
  
}

#' generate_feature_from_MS2 is to assemble MS2's precursor into features
#' 
#' @param db.MS2 results from function align_MS2_to_MS1
#' @param features results from getPeaks_multipleFiles
#' @param param ParamSet object
#' @return list as features
generate_feature_from_MS2 = function(db.MS2, features, param = new("ParamSet")){
  # convert precursor to peaks
  absMz = param@absmz
  absRt = param@absRt
  idx_umapped = which(is.na(db.MS2$MS2_to_MS1[,"MS1"]))
  unmapped_MS2 = db.MS2$MS2_to_MS1[idx_umapped,]
  unmapped_MS2_peaks = cbind(mzmin=unmapped_MS2[,"precursorMZ"], mzmax=unmapped_MS2[,"precursorMZ"],
                             rtmin=unmapped_MS2[,"retentionTime"], rtmax=unmapped_MS2[,"retentionTime"])
  rownames(unmapped_MS2_peaks) = rownames(unmapped_MS2)
  # group peaks to features
  unmapped_MS2_features = getPeakRange_multipleFiles(peaks=unmapped_MS2_peaks, cores=2, param=param, prefix = "MT")
  
  # make MS2_id-Feautre_id mapping
  MS2_ids = unmapped_MS2_features$group_FT[,"id"]
  
  idMapping = lapply(1:length(MS2_ids), function(id){
    id = MS2_ids[id]
    FT_id=names(id); 
    MS2_id=unlist(strsplit(id, split=",")); 
    result = data.frame(id=MS2_id, FT_id=FT_id);
    rownames(result) = MS2_id
    result
  })
  idMapping = do.call(rbind, idMapping)
  or = order(as.numeric(idMapping[,"id"]))
  idMapping = idMapping[or,]
  
  # refresh db.MS2$MS2_to_MS1
  db.MS2$MS2_to_MS1[idx_umapped,"MS1"] = idMapping[,"FT_id"]
  
  # refresh unmapped_MS2_features$group_FT[,"id"]
  group_MT = apply(unmapped_MS2_features$group_FT, 1, function(x){
    id = strsplit(x["id"], split=",")
    id = as.numeric(unlist(id))
    Mid = rownames(unmapped_MS2)[id]
    Mid = paste(Mid, collapse = ",")
    result = data.frame(mz=x["mz"], mzmin=x["mzmin"], mzmax=x["mzmax"], rt=x["rt"], rtmin=x["rtmin"], rtmax=x["rtmax"],
                        id=Mid, npeaks=x["npeaks"], source="MS2")
    return(result)
    
  })
  group_MT = do.call(rbind, group_MT)
  
  group_FT = features$group_FT[,c("mz","mzmin", "mzmax", "rt", "rtmin", "rtmax", "id", "npeaks")]
  group_FT = cbind(group_FT, source="MS1")
  class(group_FT) = "matrix"
  
  result = list(db.MS2=db.MS2, features=list(peaks=features$peaks, group_FT = rbind(group_FT, group_MT)))
  return(result)
}

# ----- For single Feature ----- #
#' Extract peaks mapping to certain feature
#' 
#' @param idx index of features
#' @param features list results from function getPeaks_multipleFiles, containing two data frame: peaks and features
#' @return data frame, peaks corresponding idx feature
#' @export
checkFT = function(features, idx){
  sel = features$group_FT[idx,,drop=FALSE]
  ipeaks = strsplit(sel[,"id"], split=",")
  ipeaks = as.numeric(unlist(ipeaks))
  
  peaks = features$peaks[ipeaks,,drop=FALSE]
  print(sprintf("there are %s peaks mapping to %sth feature", nrow(peaks), idx))
  return(peaks)
}

#' Compute similarity score of two sets of MS2 spectrum
#' 
#' @importFrom parallel mclapply makeCluster parLapply detectCores
#' @param MS2_set1 list containing multiple MS2 spectrum
#' @param MS2_set2 list containing multiple MS2 spectrum
#' @param cores threads for parallel computing
#' @return similarity matrix
#' @export
get_MS2_cor = function(MS2_set1, MS2_set2, cores=2){
  
  # score = mclapply(MS2_set1, function(x){
  #   spec.exp = as.matrix(x); 
  #   colnames(spec.exp)[1:2] = c("mz","intensity")
  #   
  #   # compare each with MS2 spectrum in MS2_set2 with spec.exp
  #   score1 = lapply(MS2_set2, function(y){
  #     if (is.null(y)) return(0)
  #     spec.lib = as.matrix(y)[,1:2]; 
  #     colnames(spec.lib)[1:2] = c("mz","intensity")
  #     if (nrow(spec.exp)==0 | nrow(spec.lib)==0) return(0)
  #     score = GetMatchResult(spec.exp=spec.exp, spec.lib=spec.lib)[[1]]
  #     return(score)
  #   })
  #   score2 = unlist(score1)
  #   return(score2)
  # }, mc.cores=cores)
  
  if (cores >= detectCores()) cores=detectCores()-1
  cl = makeCluster(cores)
  
  score = parLapply(cl, MS2_set1, function(x){
    spec.exp = as.matrix(x);
    colnames(spec.exp)[1:2] = c("mz","intensity")

    # compare each with MS2 spectrum in MS2_set2 with spec.exp
    score1 = lapply(MS2_set2, function(y){
      if (is.null(y)) return(0)
      spec.lib = as.matrix(y)[,1:2];
      colnames(spec.lib)[1:2] = c("mz","intensity")
      if (nrow(spec.exp)==0 | nrow(spec.lib)==0) return(0)
      score = GetMatchResult(spec.exp=spec.exp, spec.lib=spec.lib)[[1]]
      return(score)
    })
    score2 = unlist(score1)
    return(score2)
  })
  
  score = do.call(rbind, score)
  return(score)
}

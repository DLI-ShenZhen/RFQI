#' Extract spectrum from mzXML files
#' 
#' getPeaks_L is to get peaks and/or MS2 based on peak_range
#' if have ref_file, peaks and MS2 is RT-corrected
#' 
#' @param file path of one mzxml file
#' @param peak_range dataframe, each row reprensents a peak/feature area, commonly result from \code{getPeakRange_multiples}
#' @param ref_file path of referece mzxml file, using to do retention time alignment
#' @param MS1 logical, if TRUE, extract TIC intensity of peaks in peak_range
#' @param step numeric, use for integration
#' @param MS2 logical, if TRUE, extract MS2 spectrum
#' @return list
#' @export
#' @examples
#' \donttest{
#' files = files <- dir(system.file(package="RFQI", dir="extdata"), full.name=TRUE, pattern="mzXML$")
#' data("features", package="RFQI")
#' peak_range = features$group_FT
#' getPeaks_L(file=files[1], ref_file=files[2], peak_range=peak_range, MS1=TRUE, step=0.1, MS2=FALSE)
#' } 
#'  
getPeaks_L = function(file, peak_range=NULL,
                      # time alignment
                      ref_file = NULL,
                      # get MS1 peaks
                      MS1 = TRUE,
                      step = 0.1,
                      # get MS2
                      MS2 = FALSE
){
  # get adjustRtime, result is file.rtcor
  # if (!is.na(ref_file) & md5sum(file)!=md5sum(ref_file)){
  use.refFile = TRUE
  if (is.null(ref_file)) use.refFile=FALSE
  if (!is.null(ref_file)){
    if (md5sum(file)==md5sum(ref_file)) use.refFile=FALSE
  }
  
  if (!use.refFile){
    f.in = file
    fileIdx = 1
  } else {
    f.in = c(ref_file, file)
    fileIdx = 2
  }
  
  raw_data = readMSData(files=f.in, mode="onDisk", msLevel. = NULL, centroided.=TRUE)
  
  if (use.refFile){
    rtcor = adjustRtime(raw_data, param=ObiwarpParam(binSize=0.25, centerSample = 1, localAlignment = FALSE), msLevel=1L)
    rtcor = cbind(featureData(raw_data)@data, rtcor=rtcor)
  } else {
    rtcor = cbind(featureData(raw_data)@data, rtcor=featureData(raw_data)@data[,"retentionTime"])
  }
  
  file.rtcor = rtcor[rtcor$fileIdx == fileIdx,,drop=FALSE]
  file.rtcor[,"retentionTime"] = file.rtcor[,"rtcor"]
  file.rtcor = file.rtcor[,-ncol(file.rtcor)]
  file.rtcor[,"fileIdx"] = 1
  
  rm(f.in, raw_data)
  
  # --- get peaks according to peak-range ---- #
  if (MS1){
    peak_range = peak_range[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE]
    peak_range = as.matrix(peak_range)
    mode(peak_range) = "numeric"
    
    file.xcmsRaw = xcmsRaw(filename = file)
    file.xcmsRaw@scantime = file.rtcor[file.rtcor$msLevel==1,"retentionTime", drop=TRUE] 
    
    peaks = getPeaks(file.xcmsRaw, peakrange=peak_range, step=step)
    peaks[,"mz"] = apply(peaks[,c("mzmin","mzmax"),drop=FALSE], 1, median)
    peaks[,"rt"] = apply(peaks[,c("rtmin","rtmax"),drop=FALSE], 1, median)
    peaks = cbind(peaks, sample=1)
    # peaks = subset(peaks, subset=peaks[,"into"]!=0)
  } else {
    peaks = NA
  }
  
  # --- get MS2 ----- #
  if (MS2){
    feature_table_MS2 = file.rtcor[file.rtcor$msLevel==2,,drop=FALSE]
    aa = mzR::openMSfile(file)
    allPeaks_MS2 = mzR::peaks(aa)[feature_table_MS2[,"spIdx"]]
  } else {
    allPeaks_MS2 = NA
    feature_table_MS2 = NA
  }
  result = list(peaks=peaks, allPeaks_MS2=allPeaks_MS2, feature_table_MS2=feature_table_MS2, timeline=rtcor[,c("fileIdx","retentionTime","rtcor")])
  return (result)
}

#' getMS2FromFiles
#' 
#' get MS2 from multiple mzXML files
#' 
#' @param files path of mzxml files
#' @param ref_file path of referece mzxml file, using to do retention time alignment
#' @return list, MS2 is a list containing MS2 spectrum; precursorInfo is a dataframe, the ith row of precursorInfo is the ith MS2's precursor information
#' @export 
getMS2FromFiles = function(files, ref_file=NULL){
  
  MS2 = list()
  for (file in files){
    MS2[[file]] = getPeaks_L(file=file, ref_file=ref_file, MS1=FALSE, MS2=TRUE)
  }
  
  allPeaks_MS2 = c()
  for (i in MS2){
    allPeaks_MS2 = c(allPeaks_MS2, i$allPeaks_MS2)
  }
  
  feature_table_MS2 = lapply(1:length(MS2), function(i){
    # sampleName = basename(names(MS2))[i]
    sampleName = names(MS2)[i]
    df = MS2[[i]]$feature_table_MS2
    df = df[,c("precursorMZ", "retentionTime", "basePeakMZ", "basePeakIntensity", "spIdx")]
    df = cbind(df, sample=sampleName)
    return(df)
  })
  feature_table_MS2 = do.call(rbind, feature_table_MS2)
  
  return(list(MS2=allPeaks_MS2, precursorInfo=feature_table_MS2))
}

#' alignMS2_to_MS1
#' 
#' alignMS2_to_MS2 is to align MS2 database to MS1 index according to ["mzmin","mzmax","rtmin","rtmax"]
#' 
#' @param ms2Info result of getMS2FromFiles
#' @param features result of getPeakRange_multipleFiles
#' @param cores numeric, for parallel computing
#' @return list
#' @export
align_MS2_to_MS1 = function(ms2Info, features, cores=3){
  if (cores >= availableCores()) cores=availableCores()-1
  cl = makeCluster(cores)
  
  options(stringsAsFactors = FALSE)
  # peak_range = features$group_FT
  peak_range = features$group_FT[,c("mzmin","mzmax","rtmin","rtmax")]
  peak_range = cbind(peak_range, id_df=1:nrow(peak_range))    # add id_df
  mode(peak_range) = "numeric"
  # sort peak_range based on mzmin
  peak_range = peak_range[order(peak_range[,"mzmin"]),]
  # make anchor point
  if (nrow(peak_range) > 100){
    splice = nrow(peak_range) %/% 100
    anchor = 1:splice*100 + 1
    if (anchor[length(anchor)] >= nrow(peak_range)){anchor=anchor[1:(length(anchor)-1)]}
    anchor_df = peak_range[anchor,,drop=FALSE]
  }
  
  feature_table = ms2Info$precursorInfo[,c("precursorMZ", "retentionTime")]
  feature_table = as.matrix(feature_table)
  mode(feature_table) = "numeric"
  # feature_table = cbind(feature_table, id=1:nrow(feature_table))
  # feature_table = as.matrix(feature_table); mode(feature_table)="numeric"
  # Match = mclapply(1:nrow(feature_table), function(idx){
  #   ref = feature_table[idx,,drop=TRUE]
  #   
  #   # according to m/z difference, locate nearest peaks
  #   anchor_filter = abs(anchor_df[,"mzmin",drop=TRUE] - ref["precursorMZ"])
  #   anchor_point = which(anchor_filter == min(anchor_filter))
  #   if (anchor_point==1){
  #     anchor_area = 1:anchor[2]
  #   } else if (anchor_point==length(anchor)){
  #     anchor_area = anchor[(anchor_point-1)]:nrow(peak_range)
  #   } else {
  #     anchor_area = anchor[(anchor_point-1)]:anchor[(anchor_point+1)]
  #   }
  #   
  #   df = peak_range[anchor_area,]
  #   
  #   mz_filter = apply(df, 1, function(x){ref["precursorMZ"] >= x["mzmin"] & ref["precursorMZ"] <= x["mzmax"]})
  #   mz_filter = df[mz_filter,,drop=FALSE]
  #   if (nrow(mz_filter)>0){
  #     rt_filter = apply(mz_filter, 1, function(x){ref["retentionTime"] >= x["rtmin"] & ref["retentionTime"]<=x["rtmax"]})
  #     # rt_filter = ifelse(sum(rt_filter)>0, mz_filter[rt_filter,"id_df",drop=TRUE], NA)
  #     rt_filter = ifelse(sum(rt_filter)>0, rownames(mz_filter)[rt_filter], NA)
  #   } else {
  #     rt_filter = NA
  #   }
  #   return(rt_filter)
  # }, mc.cores=cores)
  
  Match = parLapply(cl, 1:nrow(feature_table), function(idx){
    ref = feature_table[idx,,drop=TRUE]

    # according to m/z difference, locate nearest peaks
    anchor_filter = abs(anchor_df[,"mzmin",drop=TRUE] - ref["precursorMZ"])
    anchor_point = which(anchor_filter == min(anchor_filter))
    if (anchor_point==1){
      anchor_area = 1:anchor[2]
    } else if (anchor_point==length(anchor)){
      anchor_area = anchor[(anchor_point-1)]:nrow(peak_range)
    } else {
      anchor_area = anchor[(anchor_point-1)]:anchor[(anchor_point+1)]
    }

    df = peak_range[anchor_area,]

    mz_filter = apply(df, 1, function(x){ref["precursorMZ"] >= x["mzmin"] & ref["precursorMZ"] <= x["mzmax"]})
    mz_filter = df[mz_filter,,drop=FALSE]
    if (nrow(mz_filter)>0){
      rt_filter = apply(mz_filter, 1, function(x){ref["retentionTime"] >= x["rtmin"] & ref["retentionTime"]<=x["rtmax"]})
      # rt_filter = ifelse(sum(rt_filter)>0, mz_filter[rt_filter,"id_df",drop=TRUE], NA)
      rt_filter = ifelse(sum(rt_filter)>0, rownames(mz_filter)[rt_filter], NA)
    } else {
      rt_filter = NA
    }
    return(rt_filter)
  })
  
  Match = unlist(Match)
  # save(Match, file="Match.rData")
  
  result = cbind(ms2Info$precursorInfo, MS1=Match)
  return(list(MS2=ms2Info$MS2, MS2_to_MS1=result))
  
}





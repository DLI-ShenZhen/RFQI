#' ID conversion, used to encoding rownames and colnames
#'
#' convert a number into a character with a prefix in fixed length
#'
#' @param x numeric
#' @param prefix character, this will be added to the head of x
#' @param n the result length of character
#' @return a character whose length is \code{n}, and with a \code{prefix}
#' @export
#' @examples
#' convert_num_to_char(1, prefix="FT", n=4)
convert_num_to_char = function(x, prefix="FT", n=4){
  x = as.character(x)
  len = nchar(x)
  imp = rep("0", n-len)
  imp = paste(imp, collapse="")
  char = paste(prefix,imp, x, sep="",collapse="")
  return(char)
}

#' ID conversion
#'
#' convert a character to a numeric
#' @param x character
#' @param prefix prefix of x
#' @return a numeric
#' @export
#' @examples
#' convert_char_to_num("FT0001", prefix="FT")
convert_char_to_num = function(x, prefix="FT"){
  num = unlist(strsplit(x, split=prefix))[2]
  num = as.numeric(num)
  return(num)
}

# ---- detect peaks ---- #
## Paramset
#' @aliases ParamSet
#' @title parameter set for peak detection and peaks group
#' 
#' @description Objects of the type \code{ParamSet} allow xcms to find peaks
#' @slot peakwith numeric(2): vector defining peakwith in Centwave method
#' @slot ppm numeric(1): numeric defining ppm in Centwave method, usually MassSpectrometry ppm
#' @slot noise numeric(1): default is 0
#' @slot absMz numeric(1)
#' @slot absRt numeric(1)
#' @slot binSize numeric(1): used for alignment
#' @rdname ParamSet
setClass("ParamSet", slots=list(peakwidth="numeric", ppm="numeric", noise="numeric", absMz="numeric", absRt="numeric", binSize="numeric"),
         prototype = list(peakwidth=c(2,30), ppm=10, noise=0, absMz = 0.015, absRt=50, binSize=0.25))


#' get peaks from multiple samples(f.in)
#'
#' detect peaks from multiple sample by centwave method
#' @importFrom xcms CentWaveParam ObiwarpParam
#' @importFrom MSnbase readMSData selectFeatureData centroided
#' @importFrom xcms findChromPeaks
#' @importFrom xcms chromPeaks adjustRtime
#' @param f.in a vertor contains paths of mzxml files
#' @param ref Path of reference sample for rentention time alignment
#' @param param ParamSet class, must contains ppm, peakwidth, and noise, these three parameter are for peak detection Using Centwave method
#' @return  if has reference sample(ref), return RT-correction peaks(["mzmin","mzmax","rtmin","rtmax"])
#' @export
#' @examples 
#' files <- dir(system.file(package="RFQI", dir="extdata"), 
#'              full.name=TRUE, 
#'              pattern="mzXML$")
#' param = new("ParamSet", binSize=0.25)
#' peaks = findPeaks_L(f.in=files, ref=NULL, param=param)
findPeaks_L = function(f.in, ref = NULL, param = new("ParamSet")){

  cwp = CentWaveParam(peakwidth=param@peakwidth, ppm=param@ppm, noise=param@noise)   ## Using CentWave algorithm to find peaks
  obi = ObiwarpParam(binSize=param@binSize, centerSample = 1, localAlignment = FALSE)
  
  if (!is.null(ref)){f.in = c(ref, f.in)}

  # Read data and find peaks using centwave
  raw_data = readMSData(files=f.in, mode="onDisk", msLevel. = NULL, centroided.=TRUE)
  xdata = findChromPeaks(raw_data, param=cwp)

  # if has reference sample, do retention time correction
  if (is.null(ref)){
    return(chromPeaks(xdata))
  } else {
    # Alignment retention time using obiwarp algorithm
    binSize = 0.6
    xdata = adjustRtime(xdata, param=obi, msLevel=1L)
    peaks_adjustRT = chromPeaks(xdata)
    peaks_adjustRT = peaks_adjustRT[!duplicated(peaks_adjustRT),]
    return(peaks_adjustRT)
  }

}

#' Whether peak area overlap
#'
#' Function judegPeakOverlap is used to judge whether two peaks ([mzmin, mzmax, rtmin, rtmax]) are overlap
#' @aliases judgePeakOverlap
#' @param ref,alt a numeric vector, orderd in "mzmin","mzmax","rtmin","rtmax"
#' @param Mz logical, whether judge m/z overlap
#' @param Rt logical, whether judge rt overlap
#' @param absMz tolerance of molecular weight difference in Da unit between ref m/z and alt m/z
#' @param absRt tolerance of retention time difference between ref rt and alt rt
#' @return logical, indicating whether ref and alt peaks are overlap
#' @export
#' @examples 
#' ref = c(70.001, 70.003, 26, 40)
#' alt = c(70.002, 70.004, 50,60)
#' judgePeakOverlap(ref=ref, alt=alt, Mz=TRUE, Rt=TRUE, absMz=0, absRt=15) 
#' 
judgePeakOverlap = function(ref, alt, Mz=TRUE, Rt=TRUE, absMz=0, absRt=15){
  # if Mz, exec Mz_filter; if Rt, exec Rt_filter
  names(ref) = names(alt) = c("mzmin", "mzmax", "rtmin", "rtmax")
  mode(ref) = mode(alt) = "numeric"

  # ---- Rt filter ---- #
  Rt_filter = TRUE
  if (Rt){
    or = order(c(ref["rtmin"], alt["rtmin"])) # sort ref and alt based on rtmin
    if (or[1]==1){r1 = ref; r2 = alt} else {r1 = alt; r2 = ref} # r1 contains the lower rtmin, r2 the larger
    Rt_filter = abs(r2["rtmin"] - r1["rtmax"]) <= absRt | r2["rtmin"] < r1["rtmax"]
  }

  # ---- Mz filter ----- #
  Mz_filter = TRUE
  if (Mz){
    or = order(c(ref["mzmin"], alt["mzmin"])) # sort ref and alt based on rtmin
    if (or[1]==1){r1 = ref; r2 = alt} else {r1 = alt; r2 = ref} # r1 contains the lower rtmin, r2 the larger
    Mz_filter = abs(r2["mzmin"] - r1["mzmax"]) <= absMz | r2["mzmin"] < r1["mzmax"]
  }

  return(Mz_filter & Rt_filter)
}

#' find overlap peaks in a mzmin-sorted data frame
#' 
#' group_FT.FUN is to compare ith peaks with peaks after it in a peakTable to find overlap peaks
#' 
#' @param i numeric, row index, compare ith peak with (i+)th peaks
#' @param df data.frame or matrix, each row is a peak area, colnames are c("mzmin","mzmax","rtmin","rtmax"), element must be numeric, df is sorted by "mzmin"
#' @param absMz,absRt used for \code{\link{judgePeakOverlap}}
#' @return a logical vector
#' @export
group_FT.FUN = function(i, df, absMz=0.01, absRt=15){
  # df is sort by mzmin
  ## group peaks which have overlap "mzmin"-"mzmax" and "rtmin"-"rtmax"
  df = cbind(df, id_df = 1:nrow(df))
  ref = df[i,c("mzmin","mzmax","rtmin","rtmax"),drop=TRUE]; mode(ref) = "numeric"
  if (i==nrow(df)){
    return (c(i))
  } else {
    df = df[(i+1):nrow(df),,drop=FALSE]
  }

  # find remote match peaks
  if (nrow(df) > 20){
    portion = ceiling(sqrt(nrow(df)))
    splice = nrow(df) %/% portion
    anchor = 1:splice*portion + 1
    if (anchor[length(anchor)] >= nrow(df)){anchor=anchor[1:(length(anchor)-1)]}
    anchor_df = df[anchor,,drop=FALSE]

    # anchor filter
    anchor_filter = apply(anchor_df, 1, function(x){
      mode(x) = "numeric"
      result = ifelse(x["mzmin"]>ref["mzmin"], x["mzmin"]-ref["mzmax"], ref["mzmin"]-x["mzmax"])
      return(result)
    })

    anchor_point = which(anchor_filter < absMz)
    if (length(anchor_point)>0){anchor_point = anchor_point[length(anchor_point)]}

    if (length(anchor_point)==0){
      anchor_area = 1:anchor[1]
    } else if (anchor_point == length(anchor)){
      anchor_area = 1:nrow(df)
    } else {
      anchor_point = anchor_point[length(anchor_point)]
      anchor_area = 1:anchor[anchor_point+1]
    }

    df = df[anchor_area,]
  }

  # compare ith peaks and (i+1):nrow(df) peaks
  mz_flag = sapply(1:nrow(df), function(j){
    alt=df[j,c("mzmin","mzmax","rtmin","rtmax"),drop=TRUE]
    mode(alt) = "numeric"
    judgePeakOverlap(ref=ref,
                     alt=alt,
                     Mz=TRUE, Rt=FALSE, absMz=absMz, absRt=absRt)
  })
  # mz_flag = sapply((i+1):nrow(df), function(j){abs(ref["mz"] - as.numeric(df[j,"mz"])) <= absMz} )
  # mz_flag = sapply((i+1):nrow(df), function(j){
  #   alt = df[j,c("mzmin","mzmax")]; mode(alt)="numeric"
  #   or = order(c(ref["mzmin"], alt["mzmin"]))
  #   if (or[1]==1){r1 = ref; r2 = alt} else {r1 = alt; r2 = ref}
  #   mz_filter = abs(r2["mzmin"] - r1["mzmax"]) <= absMz | r2["mzmin"] < r1["mzmax"]
  #   # mz_filter = r2["mzmin"] <= r1["mzmax"]
  #   return(mz_filter)
  # })
  if (length(which(mz_flag))==0) return(c(i))

  mz_filter_data = df[1:nrow(df),,drop=FALSE][mz_flag,,drop=FALSE]
  rt_flag = sapply(1:nrow(mz_filter_data), function(j){
    alt=mz_filter_data[j,c("mzmin","mzmax","rtmin","rtmax"),drop=TRUE]
    mode(alt) = "numeric"
    judgePeakOverlap(ref=ref,
                     alt=alt,
                     Mz=FALSE, Rt=TRUE, absMz=absMz, absRt=absRt)
  })
  temp_group = c(i, as.numeric(mz_filter_data[rt_flag,"id_df",drop=TRUE]))
  return(temp_group)
}

#' group overlap peaks in a peak table
#' getPeakRange_singleFiles is to group peaks from a peaks table
#' 
#' @importFrom parallel mclapply makeCluster parLapply detectCores
#' @param peaks matrix
#' @param cores numeric, for parallel computing
#' @param param ParamSet object
#' @return list
#' @export
getPeakRange_singFile = function(peaks, cores=2, param=new("ParamSet")){
  # peaks = peaks_adjustRT[1:3000,]
  absMz = param@absMz
  absRt = param@absRt
  
  if (cores >= detectCores()) cores=detectCores()-1
  cl = makeCluster(cores)
  
  if (!("id" %in% colnames(peaks))){peaks = cbind(peaks, id=1:nrow(peaks))}

  or = order(as.numeric(peaks[,"mzmin"]))
  peaks = peaks[or,]

  # group_FT = mclapply(1:nrow(peaks), group_FT.FUN, df=peaks, absMz=absMz, absRt=absRt, mc.cores=cores)
  group_FT = parLapply(cl, 1:nrow(peaks), group_FT.FUN, df=peaks, absMz=absMz, absRt=absRt)
  
  temp = group_FT
  unique_index = which(sapply(group_FT, length)==1)
  unique = temp[unique_index]
  multiple = temp[-unique_index]

  unique_FT = unique[!(unlist(unique) %in% unlist(multiple))]
  multiple_FT = list()

  while(length(multiple) != 0){
    ref = multiple[[1]]
    combined_idx = c(1)

    increment = TRUE
    while (increment){
      flag = sapply(multiple, function(alt){
        length(intersect(ref, alt)) > 0
      })
      flag = setdiff(which(flag), combined_idx)
      if (length(flag)==0){
        increment = FALSE
      } else {
        combined_idx = c(combined_idx, flag)
        ref = unique(unlist(multiple[combined_idx]))
      }
    }
    multiple = multiple[-combined_idx]
    multiple_FT[[length(multiple_FT)+1]] = unlist(ref)
  }
  group_FT = c(unique_FT, multiple_FT)

  combine_FT = lapply(group_FT, function(index, df=peaks){
    subset = df[index,,drop=FALSE]
    FT = subset[1,,drop=TRUE]
    FT["mzmin"] = min(as.numeric(subset[,"mzmin"]))
    FT["mzmax"] = max(as.numeric(subset[,"mzmax"]))
    FT["rtmin"] = min(as.numeric(subset[,"rtmin"]))
    FT["rtmax"] = max(as.numeric(subset[,"rtmax"]))
    FT["mz"] = mean(c(as.numeric(FT["mzmin"]), as.numeric(FT["mzmax"])))
    FT["rt"] = mean(c(as.numeric(FT["rtmin"]), as.numeric(FT["rtmax"])))
    FT["id"] = paste0(subset[,"id"], collapse = ",")
    FT["npeaks"] = length(unlist(strsplit(FT["id"], split=",")))

    return(FT)
  })
  combine_FT = do.call(rbind, combine_FT)

  # sort combine_FT based on m/z
  or = order(as.numeric(combine_FT[,"mzmin"]))
  combine_FT = combine_FT[or,]

  return(combine_FT)
}

#' group peaks from multiple files 
#' 
#' getPeakRange_multipleFiles is to group peaks from multiple files or a given peaks table
#' 
#' @param f.in input
#' @param peaks data frame
#' @param ref character
#' @param param ParamSet object
#' @param cores mc.cores
#' @param prefix prefix
#' @return feature table
#' @export
getPeakRange_multipleFiles = function(f.in=NULL, peaks=NULL, ref=NULL, param = new("ParamSet"), cores=2, prefix="FT"){

  print("beginning to get peaks from each sample")
  # get peak-range. if peaks is not null, return peaks; else extract peaks from multiple files
  absMz = param@absMz
  absRt = param@absRt
  
  if(is.null(peaks)){
    peaks_adjustRT = findPeaks_L(f.in = f.in, ref=ref,param=param)
  } else{
    peaks_adjustRT = peaks
  }

  print("beginning to group/combine peaks")
  # ------ combined peaks found using centwave accoding to mz and rt overlap, result is peak_range --- #
  if (!("id" %in% colnames(peaks_adjustRT))) {
    peaks_origin = cbind(peaks_adjustRT, id=1:nrow(peaks_adjustRT))     # add an id for each peak
  } else {
    peaks_origin = peaks_adjustRT
  }
  or = order(as.numeric(peaks_origin[,"mzmin"]))   # sort peaks based on mzmin. Group peaks maily based on m/z.
  peaks_origin_sort = peaks_origin[or,]            # This step is to ensure peaks only compare nearby peaks.


  # split data for parallel computation
  # small len will cost less time, e.g. when 1000 items split 4 part (len=250), group_FT cost 10s
  # while len=2, 15s; len=1000, 28s
  len = 250 # combine 14,000 peaks cost 16 minutes
  sections = ceiling(nrow(peaks_origin_sort)/len)
  f = rep(seq(1:sections), each=len)
  f = f[1:nrow(peaks_origin_sort)]
  split_data = split.data.frame(peaks_origin_sort, f=f)

  # group peaks of each split, return combine_FT
  startTime = Sys.time()
  group_FT.split_data = lapply(split_data, getPeakRange_singFile, cores=cores, param=param)
  endTime = Sys.time()
  duration = (as.numeric(endTime) - as.numeric(startTime))/60
  print(sprintf("group %s peaks divided into %s part cost %s minutes", nrow(peaks_origin_sort), length(split_data), duration))

  combine_FT = do.call(rbind, group_FT.split_data)
  or = order(as.numeric(combine_FT[,"mzmin"]))
  combine_FT = combine_FT[or,]

  # group peaks of combine_FT, return group_FT
  startTime = Sys.time()
  # the comment block can also combine/remove peaks, but consume more time
  group_FT = getPeakRange_singFile(peaks=combine_FT, cores=cores, param=param)
  # flag = TRUE
  # while (flag){
  #   rectUni = group_FT[,c("mzmin","mzmax","rtmin","rtmax")]; mode(rectUni)="numeric"
  #   rectUni = rectUnique(rectUni);
  #   flag = length(which(!rectUni)) > 0
  #   group_FT = getPeakRange_singFile(peaks=group_FT, cores=cores, absMz=absMz, absRt=absRt)
  #   or = order(as.numeric(group_FT[,"mzmin"]))
  #   group_FT = group_FT[or,]
  #
  # }

  # remove/combine overlap features
  # group_FT = combine_FT
  step = 4
  n_anchor_point = nrow(group_FT)
  flag = TRUE
  max_iter = 100
  while (flag){
    rectUni = group_FT[,c("mzmin","mzmax","rtmin","rtmax")]; mode(rectUni)="numeric"
    rectUni = rectUnique(rectUni);

    anchor_point = which(!rectUni)
    if (length(anchor_point) < n_anchor_point){
      n_anchor_point=length(anchor_point)
    } else {
      step = step + 1
    }

    if (length(anchor_point)==0){
      flag=FALSE
    } else {
      # in order to improve speed, only select overlap peaks to combine
      anchor_point = lapply(anchor_point, function(x){(x-step):(x+step)});
      anchor_point = unlist(anchor_point)
      anchor_point_filter = lapply(anchor_point, function(x){x>=1 | x<=nrow(group_FT)})
      anchor_point = anchor_point[unlist(anchor_point_filter)]
      anchor_point = unique(anchor_point)

      unique_group = group_FT[-anchor_point,,drop=FALSE]
      overlap_group = group_FT[anchor_point,,drop=FALSE]

      overlap_group = getPeakRange_singFile(peaks=overlap_group, cores=cores, param=new("ParamSet", absMz=0, absRt=0))
      group_FT = rbind(unique_group, overlap_group)

      or = order(as.numeric(group_FT[,"mzmin"]))
      group_FT = group_FT[or,]
    }

    max_iter = max_iter - 1
    if (max_iter==0) flag=FALSE

  }

  endTime = Sys.time()
  duration = (as.numeric(endTime) - as.numeric(startTime))/60
  print(sprintf("group %s peaks divided into 1 part cost %s minutes", nrow(combine_FT), duration))

  rownames(group_FT) = sapply(1:nrow(group_FT), convert_num_to_char, n=5, prefix=prefix)

  return(list(group_FT=group_FT, peaks=peaks_origin))

}

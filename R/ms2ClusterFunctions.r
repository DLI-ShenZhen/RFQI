#' Get feature index that has MS2
#' 
#' @param MS2DB db.MS2 object
#' @param n minumum number of feature should have
#' @return feature index that has MS2
#' @export
#' @examples 
#' \donttest{
#'   data("db.MS2", package="RFQI")
#'   getMS2Feature(db.MS2)
#' }
getFeatureHasMS2 = function(MS2DB, n=2){
  uniqueFeature = table(MS2DB$MS2_to_MS1[,"MS1"])
  flag = which(uniqueFeature >= n)
  filterFeature = uniqueFeature[flag]
  return(names(filterFeature))
}

#' Extract MS2 mapping to certain Feature
#' 
#' @param idx index of features
#' @param MS2DB db.MS2 object
#' @return list, similar to db.MS2, but only contains MS2 and precursorInfo of corresponding idx
#' @export
checkMS2 = function(idx, MS2DB){
  # MS2DB is db.MS2
  db.MS2 = MS2DB
  
  flag = which(db.MS2$MS2_to_MS1[,"MS1"] == idx)
  
  if (length(flag) == 0){
    warning(sprintf("there's no MS2 Spectra mapping to %sth feature", idx))
  } else {
    precursor = db.MS2$MS2_to_MS1[flag, ]
    MS2 = db.MS2$MS2[flag]
    # print(sprintf("there are %s MS2 Spectra mapping to %sth feature", length(MS2), idx))
    
    return(list(MS2=MS2, precursorInfo=precursor))
  }
}

#' Compute MS2 correlation belonging to identical feature
#' 
#' @details Compute MS2 correlation beloning to the same feature. If not appoint feature index, 
#' this function will analysis all features having MS2 in the db.MS2 object. This step is the most time-comsuming step.
#' For one feature, time complexity is O(n2), where n is the number of MS2 this feature has. We can analysis m features simultaneously,
#' and for each feature, we can compare n MS2 simultaneously. m and n corresponses to cores1 and cores2 param.
#' 
#' @param idx if assigned, only compute MS2 similarity matrix of this feature
#' @param MS2DB db.MS2 object
#' @param n minumum number of MS2 features should have, default is 2
#' @param cores1 numeric threads for idx parallel computing, e.g. simultaneously analysis cores1 feature
#' @param cores2 numeric threads for parallel computing MS2 similarity belonging to the same feature
#' @param ppm m/z tolerance of MS2
#' @param sn signal to noise ratio
#' @param maxMS2 the max number of MS2 belonging to one feature for computing similarity; if appointed, will random select maxMS2 number MS2
#' @return list, names of list are feature name, and elements are MS2 similarity matrix of MS2 belonging to the feature
#' @export
get_ms2Cor_inner = function(MS2DB, idx=NULL, n=2, cores=1, ppm=30, sn=3, maxMS2=200){
  
  if (cores >= availableCores()) cores=availableCores()-1
  cl = makeCluster(cores)
  registerPackage(cl)
  
  if (is.null(idx)) idx=getFeatureHasMS2(MS2DB=MS2DB, n=n)
  result = lapply(idx, function(id){
    MS2_set = checkMS2(id, MS2DB = MS2DB)$MS2
    if (!is.null(maxMS2)){
      if (length(MS2_set) > maxMS2) MS2_set = sample(MS2_set, maxMS2)
    }
    corMatrix = get_multiple_ms2Cor(MS2_set = MS2_set, cl=cl, ppm=ppm, sn=sn)
    rownames(corMatrix) = colnames(corMatrix) = names(MS2_set)
    corMatrix
  })
  names(result) = idx
  
  stopCluster(cl)
  
  return(result)
}


#' Select typical MS2 based on similarity matrix
#' 
#' Difference selecting typical MS2 spectrum in highest intensity or consensus MS2, here we define the typical MS2
#' is the center MS2 in similarity matrix
#' 
#' @param MS2DB db.MS2 object
#' @param corMatrix list resulting from get_ms2Cor_inner, each element is a MS2 similarity Matrix
#' @param idx if appointed, only extract typical MS2 for features in idx
#' @param include.all logical, whether extract typical MS2 for features not in corMatrix (usually contains only 1 or 2 MS2)
#' @return list, each element contains precursors' m/z and typical MS2 spectrum
#' @export
#' @examples 
#' \donttest{
#' data("db.MS2", package="RFQI")
#' corMatrix = get_ms2Cor_inner(MS2DB=db.MS2, n=3, ppm=30, sn=0, cores1=2)
#' MS2_with_mz = getTypeMS2(MS2DB = db.MS2, corMatrix = corMatrix, include.all=TRUE)
#' MS2_with_mz[1]
#' }
getTypeMS2 = function(MS2DB, corMatrix, idx=NULL, include.all = TRUE){
  
  if (!is.null(idx)) include.all=FALSE
  if (is.null(idx)) idx = names(corMatrix)
  
  typical_1 = lapply(idx, function(id){
    MS2_set = checkMS2(id, MS2DB = MS2DB)
    id.corMatrix = corMatrix[[id]]
    degree = apply(id.corMatrix, 2, sum, na.rm=TRUE) # select the center node of similarity matrix
    sel = which.max(degree)[1]
    mz = MS2_set$precursorInfo[sel,"precursorMZ"]
    MS2 = MS2_set$MS2[[sel]]
    
    list("mz"=mz, MS2=MS2)
  })
  names(typical_1) = idx
  
  if (include.all){
    idx = setdiff(unique(MS2DB$MS2_to_MS1[,"MS1"]), c(idx, NA))
    typical_2 = lapply(idx, function(id){
      MS2_set = checkMS2(id, MS2DB = MS2DB)
      mz = MS2_set$precursorInfo[1,"precursorMZ"]
      MS2 = MS2_set$MS2[[1]]
      list("mz"=mz, MS2=MS2)
    })
    names(typical_2) = typical_2
    
    typical_1 = c(typical_1, typical_2)
  }
  
  return(typical_1)
}


#' find similarity feature based on MS2 similarity
#' 
#' Compare MS2 similarity with a given m/z tolerance
#' \code{get_MS2_cor_with_mz()} compute MS2 similarity, the MS2 is filtered by m/z tolerance
#' 
#' @param MS2Set  List, each element is a list, containing MS2's precursor's m/z and MS2 spectrum, can get through \code{\link{getTypeMS2}}
#' 
#' @param cores numeric, default is 2, for parallel compute
#' @param absMz The m/z tolerance used for filtering. If two MS2's precursor m/z are within
#' tolerance, they are compared.
#' @return A matrix with m/z in the first column and separate columns for intensities
#' in the respective spectra. If peaks were merged, their m/z corresponds to the mean
#' of the two original m/z
#' @export
#' @keywords internal
get_MS2_cor_with_mz = function(MS2_set, cores=2, absMz=0.04){
  
  if (cores >= availableCores()) cores=availableCores()-1
  cl = makeCluster(cores)
  registerPackage(cl)
  
  registerParentVars(cl)
  score = parLapply(cl, MS2_set, function(x){
    mz.exp = x$mz
    spec.exp = x$MS2; colnames(spec.exp)[1:2] = c("mz","intensity")
    
    score1 = lapply(MS2_set, function(y){
      mz.lib = y$mz
      if (abs(mz.lib - mz.exp) > absMz) return(NA)
      
      spec.lib = y$MS2;
      if (is.null(spec.lib)) return(0)
      colnames(spec.lib)[1:2] = c("mz","intensity")
      if (nrow(spec.exp)==0 | nrow(spec.lib)==0) return(0)
      score1 = GetMatchResult(spec.exp=spec.exp, spec.lib=spec.lib)
      return(score1)
    })
    
    score2 = unlist(score1)
    return(score2)
  })  
  stopCluster(cl)
  
  score = do.call(rbind, score)
  diag(score) = NA
  colnames(score) = rownames(score) = names(MS2_set)
  
  return(score)
}

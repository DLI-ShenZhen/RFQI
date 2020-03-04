#' add adduct m/z to metabolite
#' 
#' @param lib spectrum library, A list
#' @param type column and ionization mode, could be one of c("RP_pos","RP_neg","HILIC_pos", "HILIC_neg")
#' @return a matrix, represents each metablites adduct m/z
#' @examples
#' \donttest{
#'   data("spectrumDB", package="RFQI")
#'   data("adductDB", package="RFQI")
#'   add_adduct(lib=spectrumDB, type="RP_pos")
#' }
add_adduct = function(lib, type="RP_pos"){
  
  data("adductDB", package = "RFQI")
  adduct = match.arg(type, c("RP_pos","RP_neg","HILIC_pos", "HILIC_neg"))
  adduct = adductDB[[adduct]]
  
  mzs = as.numeric(lib$Info[,"mz",drop=TRUE])
  lib.mz <- t(sapply(mzs, function(mz) {
    mz <- as.numeric(mz)
    apply(adduct, 1, function(info.adduct) {
      x <- gsub('\\(', '',
                gsub('M.*', '', info.adduct['adduct']))
      xm <- ifelse(x == '', 1, as.numeric(x))
      xm * mz + as.numeric(info.adduct['mz'])
    })
  }))
  rownames(lib.mz) = rownames(lib$Info)
  colnames(lib.mz) = adduct[,"adduct"]
  
  return(lib.mz)
}

#' Compute FDR
#' 
#' Compute possibility of one features generates certain (often) MS2, if possibility greater than a threshold, we say this feature is the metabolite
#' 
#' @param score similarity vector of one spectrum with multiple MS2
#' @param ms2_inner_score similarity matrix of multiple MS2
#' @param ratio ratio threshold, is the minimum possibility of feature generating certain MS2
#' @return ratio
getFDR = function(score, ms2_inner_score, ratio=0.1){
  # "every" compare each inner score with spec.lib, and using ratio(freq) agree with spec.lib to determine
  
  # ensure score is a matrix which each column is a metabolite
  if (nrow(score)==1) score = t(score)
  if (is.null(nrow(score))){score=as.matrix(score)}
  
  # compute pros/cons ratio
  FDR_filter = apply(score, 2, function(x){
    # alternative="less" is the alternative that x has a lower mean than y
    # i.e. if p.value<0.05, receive alternative, namely x < y; if p.value>0.05, x>=y
    summary_x = c(median(x,na.rm=T), max(x,na.rm=T), sd(x,na.rm=T))
    
    p.value.greater = apply(ms2_inner_score, 1, function(y){
      t.test(x, y, alternative = "less")$p.value
    })
    pros = sum(p.value.greater>0.05)
    cons = sum(p.value.greater<=0.05)
    
    flag = ifelse(pros/(pros+cons) >= ratio, 1, 0)
    return(c(summary_x, flag, pros=pros, cons=cons, ratio=pros/(pros+cons)))
  })
  FDR_filter = t(FDR_filter)  
  
  # Judge whether ratio meets cutoff
  flag = which(as.numeric(FDR_filter[,1]) == max(as.numeric(FDR_filter[,1]))) 
  
  return(FDR_filter[flag,,drop=FALSE])
  
}

#' Identify features based on multiple MS2 similarity
#' 
#' Need to prepare standard metabolites MS2 spectrum database
#' 
#' @param lib standard metabolite MS2 spectrum database, We have offer a database containing 3000 metabolites, using \code{data("spectrumDB", package="RFQI")}
#' @param MS2DB db.MS2 object, result from function align_MS2_to_MS1
#' @param MS2_inner_cor list containing ms2 similarity score matrix 
#' @param MS1_idx if not NULL, only identify appoint MS1_index, else identify all features having MS2
#' @param absMz mz tolerance in Da unit, this is used for find matched precursor
#' @param adduct possible adduct when ionization, can be one of c("RP_pos", "RP_neg", "HILIC_pos", "HILIC_neg"), RP represents reverse column, pos is positive mode
#' @param cores cores for parallel computing
#' @return identification result
#' @export
get_identify = function(lib, MS2DB, MS2_inner_cor, MS1_idx=NULL, absMz = 0.015,
                        adduct="RP_pos", cores=1){
  meta.precursor = add_adduct(lib=lib, type = adduct)
  if (cores >= availableCores()) cores=availableCores()-1
  cl = makeCluster(cores)
  
  # if appoint MS1_idx, only identify MS1_idx, else identify all features with MS2 
  if (!is.null(MS1_idx)){
    MS1_index=MS1_idx
  } else {
    MS1_index = getFeatureHasMS2(MS2DB = MS2DB, n=1) 
  }
  MS2_to_MS1 = MS2DB$MS2_to_MS1
  
  # score = mclapply(MS1_index, function(idx){
  #   
  #   MS2.idx = which(MS2_to_MS1[,"MS1"] %in% idx)
  #   MS2.idx.mz = median(MS2_to_MS1[MS2.idx,"precursorMZ"])
  #   
  #   mz.filter = abs(MS2.idx.mz - as.numeric(meta.precursor)) <= absMz
  #   mz.filter = which(mz.filter)
  #   axe.col = ceiling(mz.filter / nrow(meta.precursor))
  #   axe.row = mz.filter - (axe.col-1)*nrow(meta.precursor)
  #   
  #   lib.temp = lib$spectrum[rownames(meta.precursor)[axe.row]]
  #   
  #   idx.MS2 = MS2DB$MS2[MS2.idx]
  #   if (idx %in% names(MS2_inner_cor)) idx.MS2 = MS2DB$MS2[rownames(MS2_inner_cor[[idx]])]
  #   
  #   # compare db.MS2 and idx.MS2
  #   score = get_MS2_cor(MS2_set1 = idx.MS2, MS2_set2 = lib.temp)
  #   return(list(score=score, adduct=colnames(meta.precursor)[axe.col]))
  # }, mc.cores=cores) 
  
  score = parLapply(cl, MS1_index, function(idx){

    MS2.idx = which(MS2_to_MS1[,"MS1"] %in% idx)
    MS2.idx.mz = median(MS2_to_MS1[MS2.idx,"precursorMZ"])

    mz.filter = abs(MS2.idx.mz - as.numeric(meta.precursor)) <= absMz
    mz.filter = which(mz.filter)
    axe.col = ceiling(mz.filter / nrow(meta.precursor))
    axe.row = mz.filter - (axe.col-1)*nrow(meta.precursor)

    lib.temp = lib$spectrum[rownames(meta.precursor)[axe.row]]

    idx.MS2 = MS2DB$MS2[MS2.idx]
    if (idx %in% names(MS2_inner_cor)) idx.MS2 = MS2DB$MS2[rownames(MS2_inner_cor[[idx]])]

    # compare db.MS2 and idx.MS2
    score = get_MS2_cor(MS2_set1 = idx.MS2, MS2_set2 = lib.temp)
    return(list(score=score, adduct=colnames(meta.precursor)[axe.col]))
  })
  
  # voting2
  ID = lapply(1:length(score), function(i, cutoff=0.5){
    idx_MS1 = MS1_index[i]
    item = score[[i]]
    
    df = item$score
    # NULL filter
    if (is.null(df)) return(rep(NA,9))
    
    # cutoff filter
    flag_cutoff = apply(df, 2, function(x){mean(x>cutoff)>=0.5})  # For each column, whether half (0.5) elements greater than cutoff
    flag_cutoff.pass = sum(flag_cutoff) >= 1                           # at least one column pass cutoff filter
    if (!flag_cutoff.pass){return(rep(NA,9))}

    df = df[,which(flag_cutoff),drop=FALSE]
    # compute FDR
    ref_MS2_score = MS2_inner_cor[[idx_MS1]]

    if (is.null(ref_MS2_score)){
      flag = apply(df, 2, sum)
      flag = which(flag == max(flag))[1]
      return(c(median(df[,flag]), max(df[,flag]), sd(df[,flag]), rep(NA,4), colnames(df)[flag], item$adduct[flag]))
    } else if (nrow(ref_MS2_score)<3 | length(setdiff(unique(as.numeric(ref_MS2_score)),NA))==1){
      flag = apply(df, 2, sum)
      flag = which(flag == max(flag))[1]
      return(c(median(df[,flag]), max(df[,flag]), sd(df[,flag]), rep(NA,4), colnames(df)[flag], item$adduct[flag]))

    } else {

      FDRFilter = getFDR(score=df, ms2_inner_score = ref_MS2_score, ratio=0.1)
      flag = which(colnames(item$score)==rownames(FDRFilter))
      return(c(FDRFilter[1,,drop=TRUE], rownames(FDRFilter)[1], item$adduct[flag]))
    }
  })
  
  ID = do.call(rbind, ID); 
  colnames(ID) = c("ScoreMedian", "ScoreMax", "ScoreSD", "FDRFilter", "num_pros", "num_cons", "ratio", "labid", "adduct")
  rownames(ID) = MS1_index
  ID_filter = ID[!is.na(ID[,"ScoreMax"]),,drop=FALSE]
  
  if (nrow(ID_filter)==0){
    anno = cbind(ID_filter, metabolite=NULL)
  } else {
    # get anno of MS1
    or = order(rownames(ID_filter))
    anno = ID_filter[or,,drop=FALSE]
    
    metabolite =  lib$Info[anno[,"labid"], -1]
    anno = cbind(anno, metabolite)
  }
  names(score) = MS1_index
  
  return(list(anno=anno, score=score))
}




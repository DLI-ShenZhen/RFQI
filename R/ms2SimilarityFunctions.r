# --- MS2 similarity algorithm --- #

#' translate Da to ppm
#' 
#' Calculating difference of two m/z, then convert diff into ppm
#' 
#' @param mz numeric vector, m/z values sorted by inreasing
#' @param mz.ppm.thr numeric, threshold of ppm difference
#' @return vector, ppm difference of m/z with the first m/z
#' @examples 
#' \donttest{
#'   GetDiffMZppm(c(100.01,100.02, 100.03))
#' }
GetDiffMZppm <- function(mz, mz.ppm.thr = NULL) {
  mz.diff <- diff(mz) / mz[-1] * 1e6
  if (!is.null(mz.ppm.thr)) {
    idx <- which(mz[-1] <= mz.ppm.thr)
    mz.diff[idx] <- mz.diff[idx] * mz[-1][idx] / mz.ppm.thr
  }
  mz.diff
}


#' group m/z based on m/z tolerance
#'
#' As MassSpectrometry has fluctuate in m/z detection, one ion can be detected as two similar m/z, thus we need to combine these two m/z.
#' After Combining, m/z is the mean, and intensity is the sum.
#' @param spec numeric matrix, columns are "mz" and "intensity", represents an MS2 spectrum
#' @param ppm.ms2match numeric, m/z tolerance of MS2 in ppm mode
#' @param mz.ppm.thr used in GetDiffMZppm function
#' @return MS2 spectrum
MatchSpec <- function(spec, ppm.ms2match = 30, mz.ppm.thr = 400) {
  while (TRUE) {
    mz.diff.ppm <- GetDiffMZppm(spec[, 'mz'], mz.ppm.thr = mz.ppm.thr)
    idx <- which(mz.diff.ppm < ppm.ms2match)
    if (length(idx) > 0) {
      i <- tail(idx, 1)
      j <- which.max(spec[c(i, i + 1), 'intensity'])
      spec[i, 'intensity'] <- spec[i + j - 1, 'intensity']
      i2 <- i + 1
      spec[i, 'mz'] <- spec[i2, 'mz']
      spec <- spec[-i - 1, , drop = FALSE]
    } else {
      break
    }
  }
  return(spec)
}

#' combine two MS2 spectrum based on m/z
#' 
#' @param spec.exp MS2 spectrum
#' @param spec.lib MS2 spectrum
#' @param ppm.ms2match m/z tolerance of MS2 
#' @return list, contains m/z alignment exp and lib spectrum
GetSpec2Match <- function(spec.exp, spec.lib, ppm.ms2match = 30) {
  
  # align m/z of spec.exp and spec.lib
  mz.pool   <- sort(c(spec.exp[, 1], spec.lib[, 1]))	# combine mz
  spec.exp.pool = spec.lib.pool = cbind('mz' = mz.pool, 'intensity' = 0)		# construct alignment matrix
  spec.exp.pool[match(spec.exp[,1], spec.exp.pool), "intensity"] = spec.exp[,2]
  spec.lib.pool[match(spec.lib[,1], spec.lib.pool), "intensity"] = spec.lib[,2]
  
  # combine nearby peaks
  pk.spec  <- MatchSpec(spec.exp.pool, ppm.ms2match = ppm.ms2match)	# combine adjacent m/z if GetDiffMZppm < pps.ms2match)
  lib.spec <- MatchSpec(spec.lib.pool, ppm.ms2match = ppm.ms2match)	# combine adjacent m/z if GetDiffMZppm < pps.ms2match)
  
  return(list('exp' = pk.spec, 'lib' = lib.spec))
}


#' compute cosine similarity of two MS2 spectrum
#'
#' @param spec.exp experimental MS2 spectrum, columns are "mz" and "intensity"
#' @param spec.lib reference MS2 spectrum, columns are "mz" and "intensity"
#' @param sn signal to noise ratio
#' @param ppm m/z tolerence of MS2
#' @return cosine similarity of spec.exp and spec.lib
#' @export
GetMatchResult = function(spec.exp, spec.lib, sn=3, ppm=30){
  # spec.lib=a; spec.exp=b; sn=10
  colnames(spec.lib)=colnames(spec.exp)=c("mz","intensity")
  mode(spec.lib) = mode(spec.exp) = "numeric"
  spec.lib[,2] = spec.lib[,2]/max(spec.lib[,2])*100
  spec.exp[,2] = spec.exp[,2]/max(spec.exp[,2])*100
  spec.lib.filter = spec.lib[spec.lib[,2]>sn,,drop=FALSE]
  spec.exp.filter = spec.exp[spec.exp[,2]>sn,,drop=FALSE]
  
  # combine spec.lib and spec.exp
  spec2match <- GetSpec2Match(spec.exp, spec.lib,                               # GetSpec2Match returns match(spec.exp, spec.exp+spec.lib) for "exp"
                              ppm.ms2match = ppm)			# and match(spec.lib, spec.exp+spec.lib) + spec.exp for "lib"
  
  # compute cosine similarity
  int.spec = spec2match$exp[,"intensity"]
  int.lib = spec2match$lib[,"intensity"]
  
  score = sum(int.spec*int.lib)/sqrt(sum(int.spec^2) * sum(int.lib^2))
  
  return(score)
}

#' compute pair cosine similarity score of multiple MS2
#' 
#' @aliases get_MS2_cor_inner
#' @param MS2_set a list whose elements are MS2 spectrum
#' @param cores numeric for parallel computing
#' @param cl cluster generate by parallel::makecluster, if cl is NULL, then use cores to makeCluster
#' @param ppm m/z tolerance
#' @param sn signal to noise ratio
#' @return a similarity score matrix
#' @export
get_multiple_ms2Cor = function(MS2_set, cores=2, cl=NULL, ppm=30, sn=3){
  
  # if appoint cl, using parLapply
  # if appoint cores, using mclapply
  
  if (is.null(cl)){
    # if (cores >= availableCores()) cores=availableCores()-1
    # cl = makeCluster(cores)
    # registerPackage(cl)
    # closeCluster = TRUE
    score = mclapply(1:(length(MS2_set)-1), function(i){
      spec.exp = MS2_set[[i]]; colnames(spec.exp)[1:2] = c("mz", "intensity")
      
      # compare spec.exp with it's latter MS2
      score1 = lapply(MS2_set[(i+1):length(MS2_set)], function(y){
        if (is.null(y)) return(0)
        spec.lib = y; colnames(spec.lib)[1:2] = c("mz","intensity")
        if (nrow(spec.exp)==0 | nrow(spec.lib)==0) return(0)
        score1 = GetMatchResult(spec.exp=spec.exp, spec.lib=spec.lib, ppm=30, sn=3)
        return(score1)
      })
      score2 = unlist(score1)
      return(score2)
    }, mc.cores=cores)
    
  } else {
    registerParentVars(cl)
    score = parLapply(cl, 1:(length(MS2_set)-1), function(i){
      spec.exp = MS2_set[[i]]; colnames(spec.exp)[1:2] = c("mz", "intensity")
      
      # compare spec.exp with it's latter MS2
      score1 = lapply(MS2_set[(i+1):length(MS2_set)], function(y){
        if (is.null(y)) return(0)
        spec.lib = y; colnames(spec.lib)[1:2] = c("mz","intensity")
        if (nrow(spec.exp)==0 | nrow(spec.lib)==0) return(0)
        score1 = GetMatchResult(spec.exp=spec.exp, spec.lib=spec.lib, ppm=30, sn=3)
        return(score1)
      })
      score2 = unlist(score1)
      return(score2)
    })
  }
  
  
  result = matrix(nrow=length(MS2_set), ncol=length(MS2_set))
  for (i in 1:length(score)){
    result[i,(i+1):ncol(result)] = score[[i]]
  }
  result[lower.tri(result)] = result[upper.tri(result)]
  
  return(result)
}

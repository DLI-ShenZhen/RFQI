#' get_compareFT
#'
#' compare two feature table overlap
#' 
#' @param M1 first feature table
#' @param M2 second feature table
#' @param absMz numeric
#' @param absRt numeric
#' @param cores thread for parallel computing
#' @export
#' @return list
get_compareFT = function(M1, M2, absMz=0.005, absRt=15, cores=2){
  M1 = as.matrix(M1)[,c("mzmin","mzmax","rtmin","rtmax")]
  M2 = as.matrix(M2)[,c("mzmin", "mzmax","rtmin","rtmax")]
  mode(M1) = mode(M2) = "numeric"
  
  M2.anchor = seq(1,nrow(M2),200)
  M2.anchor[length(M2.anchor)]=nrow(M2)
  M2.anchor.df = M2[M2.anchor,,drop=FALSE]
  
  matchID = mclapply(1:nrow(M1), function(i){
    x = M1[i,,drop=TRUE]
    target1 = apply(M2.anchor.df, 1, function(y){as.numeric(x["mzmin"])-as.numeric(y["mzmin"])})
    target1 = which.min(abs(target1))
    start = ifelse(target1==1,1,target1-1)
    end = ifelse(target1==nrow(M2.anchor.df), target1, target1+1)
    
    target.df = M2[M2.anchor[start]:M2.anchor[end],,drop=FALSE]
    result = apply(target.df, 1, judgePeakOverlap, alt=x, absMz=absMz, absRt=absRt)
    names(which(result))
  }, mc.cores=cores)
  names(matchID) = rownames(M1)
  
  return(matchID)
}

#' feature unique match
#' 
#' @param M1_to_M2 results from get_compareFT, M1 feature table mapping to M2 feature table
#' @param M2_to_M1 reverse of M1_to_M2
#' @return data frame
match.one2one = function(M1_to_M2, M2_to_M1){
  obj = M1_to_M2
  
  # one M1 ID only mapping to one M2 ID
  len = sapply(obj, length); 
  obj.refresh1 = obj[len==1]
  
  # one FT.serial ID only mapping to one FT.random ID
  M2ID = unlist(obj.refresh1)
  M2ID.unique = setdiff(M2ID, M2ID[duplicated(M2ID)])
  
  obj.refresh2 = sapply(obj.refresh1, function(x){x[1] %in% M2ID.unique})
  obj.refresh2 = obj.refresh1[names(which(obj.refresh2))]
  
  # result
  result = data.frame(M1ID=names(obj.refresh2), M2ID=unlist(obj.refresh2), stringsAsFactors = FALSE)
  rownames(result) = names(obj.refresh2)
  return(result)
}

#' add error bar for boxplot
#' 
#' @param x0 numeric
#' @param y0 numeric
#' @param x1 numeric
#' @param y1 numeric
#' @param y numeric
#' @param pvalue character
#' @return null
#' @export
errorBar = function(x0=1, y0, x1=2, y1, y, pvalue){
  segments(x0 = x0, y0 = y0, x1 = x0, y1 = y)
  segments(x0 = x1, y0 = y1, x1=x1, y1=y)
  segments(x0 = x0, y0 = y, x1 = x1, y1 = y)
  text=sprintf("p < %s",pvalue)
  text(x=mean(c(x0,x1)), y=y, labels = text, adj=c(0.5,-1), cex=1.2)
  print(mean(x0,x1))
}

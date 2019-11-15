#' Do MS2 clustering using dbscan
#' 
#' @param idx index of features
#' @param MS2_inner_cor list, names are feature index, and each element is a MS2 similarity matrix
#' @param eps similarity score threshold
#' @param minPts minimum member of a node
#' @return a list, cluster contains each cluster's member, and noise are isolate MS2
#' @export
get_ms2cluster = function(idx, MS2_inner_cor, eps=0.5, minPts=2){
  # test mode
  if (FALSE){
    idx = "FT00773"
    MS2_inner_cor = MS2_inner_cor
    eps = 0.5
    minPts = 2
  }
  
  corMatrix = MS2_inner_cor[[idx]]
  # minPts = min(ceiling(0.1*nrow(corMatrix)), 2)
  
  # label corMatrix
  if (is.null(rownames(corMatrix))){rownames(corMatrix) = colnames(corMatrix) = as.character(1:nrow(corMatrix))}
  
  # --- Using dbscan to cluster --- #
  # init param
  seed_init = apply(corMatrix, 2, sum, na.rm=TRUE)
  seed_init = which(seed_init == max(seed_init))
  seed_init = rownames(corMatrix)[seed_init]
  
  kernal = c()  
  reachable = list()
  margin = c()
  init_nrow = nrow(corMatrix)
  dbscan_cluster = list()
  
  seed = seed_init
  while (TRUE){
    # minimum members of each cluster
    if (nrow(corMatrix)<=minPts){break}
    
    # kernal - reachable node
    density_reachable = which(corMatrix[,seed] >= eps)
    density_reachable = names(density_reachable)
    pts = length(density_reachable)
    if (pts>= minPts){kernal = c(kernal, seed)}  # judge whether seed is kernal node
    
    # reachable-margin node
    for (i in density_reachable){    # elements of density_reachable are rownames of corMatrix
      density_reachable_i = which(corMatrix[,i]>=eps)
      density_reachable_i = rownames(corMatrix)[density_reachable_i]
      pts_i = length(density_reachable_i)
      
      if (pts_i>=minPts){
        kernal = c(kernal, i)
      }
      
      outer = setdiff(density_reachable_i, c(density_reachable, seed))
      margin = unique(c(margin, outer))    # margin may be shared among different clusters
    }
    
    # refresh dbscan_cluster, also can save kernal and margin node list
    dbscan_cluster[[length(dbscan_cluster)+1]] = c(seed, density_reachable, unique(margin))
    
    # refresh corMatrix
    remain = setdiff(rownames(corMatrix), c(seed, density_reachable, unique(margin)))
    if (length(remain)==0) break      # there is no column unclustered
    
    filter = setdiff(rownames(corMatrix), c(seed, density_reachable))
    corMatrix = corMatrix[filter, filter, drop=FALSE]
    if (nrow(corMatrix)<init_nrow){init_nrow=nrow(corMatrix)} else {break}   # corMatrix cannot be splited
    
    # refresh seed
    seed = apply(corMatrix, 2, sum, na.rm=TRUE)
    seed = which(seed == max(seed))
    seed = rownames(corMatrix)[seed]
    
    # re-init 
    kernal = c()
    density_reachable = c()
    margin = c()
    dbscan_cluster
  }
  
  if (length(dbscan_cluster)==1) {names(dbscan_cluster)=idx} else{
    names(dbscan_cluster) = paste(idx, 1:length(dbscan_cluster), sep=".")
  }
  
  # As margin nodes are share among different clusters, we need to define which cluster they belong to
  corMatrix = MS2_inner_cor[[idx]]
  common_node = unlist(dbscan_cluster); common_node=common_node[duplicated(common_node)]
  if (length(common_node) != 0){
    # for (node in common_node){
    #   node_score = corMatrix[node,,drop=TRUE]
    #   sapply(dbscan_cluster, function(j){
    #     mean(node_score[j], na.rm=TRUE)
    #   })
    # }
    b = lapply(common_node, function(i){
      i = common_node[1]
      i_score = corMatrix[i,]
      sapply(dbscan_cluster, function(j){
        a = i_score[j]
        mean(a, na.rm=TRUE)
      })
      
    })
    b = do.call(rbind,b)
    rownames(b) = common_node
    
    b = apply(b, 1, function(x){flag=which(x==max(x)); names(x)[flag]})
    
    dbscan_cluster_rmdup = lapply(1:length(dbscan_cluster), function(i){
      x = dbscan_cluster[[i]]
      x_name = names(dbscan_cluster)[i]
      x_unique = setdiff(x, common_node)
      x_common = rownames(b)[which(b==x_name)]
      return(c(x_unique, x_common))
    })
    
    len_filter = sapply(dbscan_cluster_rmdup, length) != 0
    dbscan_cluster= dbscan_cluster_rmdup[len_filter]
    
  }
  
  # filter noise, i.e cluster which members less than 2
  noise_filter = sapply(dbscan_cluster, length) <= 2
  noise = unlist(dbscan_cluster[noise_filter])
  
  dbscan_cluster = dbscan_cluster[!noise_filter]
  dbscan_cluster
  
  if (length(dbscan_cluster)==0) {
  } else if (length(dbscan_cluster)==1) {names(dbscan_cluster)=idx} else{
    names(dbscan_cluster) = paste(idx, 1:length(dbscan_cluster), sep=".")
  }
  
  return(list(cluster=dbscan_cluster, noise=noise))
}
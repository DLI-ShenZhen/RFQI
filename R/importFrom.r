#' import functions from other packages
#' 
#' not export
#' 
#' @import xcms
#' @import pheatmap
#' @import grDevices 
#' @import graphics
#' @import RColorBrewer
#' @import future
#' @import BiocParallel
#' @importFrom MSnbase readMSData selectFeatureData centroided
#' @importFrom xcms CentWaveParam ObiwarpParam
#' @importFrom xcms findChromPeaks
#' @importFrom xcms chromPeaks adjustRtime
#' @importFrom BiocParallel register MulticoreParam SnowParam
#' @importFrom parallel mclapply makeCluster parLapply stopCluster clusterExport clusterEvalQ
#' @importFrom future availableCores
#' @importFrom scales alpha
#' @importFrom tools md5sum
#' @importFrom stats median runmed
#' @importFrom utils tail
#' @importFrom stats sd t.test
#' @importFrom utils data
#' @importFrom Biobase featureData
#' @importFrom igraph graph_from_data_frame plot.igraph
#' 
loadPackage = function(){}

platform = Sys.info()[['sysname']]



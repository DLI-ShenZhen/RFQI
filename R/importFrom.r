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
#' @importFrom MSnbase readMSData selectFeatureData centroided
#' @importFrom xcms CentWaveParam ObiwarpParam
#' @importFrom xcms findChromPeaks
#' @importFrom xcms chromPeaks adjustRtime
#' @importFrom parallel mclapply makeCluster parLapply stopCluster clusterExport clusterEvalQ
#' @importFrom future availableCores
#' @importFrom scales alpha
#' @importFrom tools md5sum
#' @importFrom stats median runmed
#' @importFrom utils tail
#' @importFrom stats sd t.test
#' @importFrom utils data
#' @importFrom Biobase featureData
#' 
loadPackage = function(){}




#' Features resulting from mzXML files in extdata dir
#' 
#' @format A list containing two data frame
#' \describe{
#'   \item{group_FT}{features after combining peaks}
#'   \item{peaks}{peaks extracted from mzXML files using CentWave method}
#' }
#' @source data("features", package="RFQI")
"features"

#' MS2 spectrum resulting from getMS2FromFiles function
#' 
#' @format A list containing two elements
#' \describe{
#'   \item{MS2}{MS2 spectrum}
#'   \item{presursorInfo}{precursor information corresponding MS2 spectrum}
#' }
#' @source data("MS2_spec", package="RFQI")
"MS2_spec"

#' MS2 spectrum Mapping to features 
#' 
#' @format A list containing two elements
#' \describe{
#'   \item{MS2}{MS2 spectrum}
#'   \item{presursorInfo}{precursor information corresponding MS2 spectrum}
#' }
#' @source data("db.MS2", package="RFQI")
"db.MS2"

#' standard metabolite MS2 spectrum
#' 
#' @format A list containing two elements
#' \describe{
#'   \item{spectrum}{standard MS2 spectrum, each spectrum is a list, whose name is the metabolite ID, and value is spectrum}
#'   \item{Info}{metabolite information, including Labid, HMDBID, mz, formula, and HMDB name}
#' }
#' @source data("spectrumDB", package="RFQI")
"spectrumDB"

#' possible adduct when ionization
#' 
#' @format A list containing four matrix, each matrix is an adduct table
#' \describe{
#'   \item{HILIC_pos}{chromatographic column is HILIC, and ionization mode is "positive"}
#'   \item{HILIC_neg}{chromatographic column is HILIC, and ionization mode is "negative"}
#'   \item{RP_pos}{chromatographic column is RP, and ionization mode is "positive"}
#'   \item{RP_neg}{chromatographic column is RP, and ionization mode is "negative"}
#' }
#' @source data("adductDB", package="RFQI") 
"adductDB"
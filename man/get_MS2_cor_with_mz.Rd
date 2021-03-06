% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms2ClusterFunctions.r
\name{get_MS2_cor_with_mz}
\alias{get_MS2_cor_with_mz}
\title{find similarity feature based on MS2 similarity}
\usage{
get_MS2_cor_with_mz(MS2_set, cores = 2, absMz = 0.04)
}
\arguments{
\item{cores}{numeric, default is 2, for parallel compute}

\item{absMz}{The m/z tolerance used for filtering. If two MS2's precursor m/z are within
tolerance, they are compared.}

\item{MS2Set}{List, each element is a list, containing MS2's precursor's m/z and MS2 spectrum, can get through \code{\link{getTypeMS2}}}
}
\value{
A matrix with m/z in the first column and separate columns for intensities
in the respective spectra. If peaks were merged, their m/z corresponds to the mean
of the two original m/z
}
\description{
Compare MS2 similarity with a given m/z tolerance
\code{get_MS2_cor_with_mz()} compute MS2 similarity, the MS2 is filtered by m/z tolerance
}
\keyword{internal}

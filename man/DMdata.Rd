% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{DMdata}
\alias{DMdata}
\title{A Dataset Simulating Single-Cell Epigenomic Data}
\format{A list with three items:
\describe{
  \item{data}{a 100 by 5000 matrix, indicating 100 observations with 5000 features. 
  Each element is a non-negative integer, which could be regarded as read counts.}
  \item{nCluster}{true number of clusters.}
  \item{true.label}{a 100-dimension integer vector indicating the cluster each observation 
  should be assigned to.}
  \item{noisy.label}{a 20-dimension integer vector indicating the positions of noisy 
  observations. 90\% of the features in the noisy observations are set to be 0, 
  simulating missing values in real experiments.}
}}
\source{
Jason D Buenrostro, \emph{et al.} (2015) 
Single-cell chromatin accessibility reveals principles of regulatory variation. 
\emph{Nature}, \bold{523(7561)}, 486-490.
}
\usage{
DMdata
}
\description{
A simulated dataset used in the `SWKM` vignette for illustration to perform weighted
K-Means clustering.
}
\examples{
data <- DMdata$data
true.label <- DMdata$true.label
noisy.label <- DMdata$noisy.label
}
\keyword{datasets}

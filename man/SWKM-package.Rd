\name{SWKM-package}
\alias{SWKM-package}
\alias{SWKM}
\docType{package}
\title{
  Sparse Weighted K-Means
}
\description{
  An implementation of (sparse) weighted K-Means clustering algorithm 
    on observations with weights.
}
\details{
  Main Functions:
  
  \code{kmeans.weight}: Perform weighted K-Means algorithm on data.

  \code{kmeans.weight.tune}: Choose weight parameter for (sparse) weighted K-Means algorithm. Usually used before \code{kmeans.weight} or \code{KMeansSparseCluster.weight}.

  \code{KMeansSparseCluster.weight}: Perform sparse weighted K-Means algorithm on data.

  \code{KMeansSparseCluster.permute.weight}: Choose sparsity parameter for sparse weighted K-Means algorithm. Usually used before \code{KMeansSparseCluster.weight}, and after weight parameter is tuned or known.

  \code{ChooseK}: Choose the number of clusters K for (sparse) weighted K-Means clustering. Usually used before clustering method is performed.
  
  Please refer to the vigenette for more details.
}
\author{
Wenyu ZHANG.

Maintainer: Wenyu ZHANG <wyzhangxii@outlook.com>
}
\references{
  Daniela M Witten and Robert Tibshirani (2010). A framework for feature selection in clustering.  \emph{Journal of the American Statistical Association}, \bold{105(490)}, 713-726.
  
  Robert, T. \emph{et al.} (2001). Estimating the number of clusters in a data set via the gap statistic. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{63(2)}, 411-423.
}
\keyword{ package, sparse, weighted, k-means }

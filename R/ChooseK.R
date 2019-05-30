#' Choose the Number of Clusters K for (Sparse) Weighted K-Means Clustering
#'
#' The number of clusters K should be determined before the clustering method is performed. 
#' A permutation approach using Gap Statistic is used.
#' @inheritParams kmeans.weight
#' @param nClusters a candidate sequence of K. Default is \code{2:6}. 
#' @param nperms Number of permutations. Default is \code{20}.
#' @param ... unused.
#' @importFrom stats kmeans
#' @keywords Weighted K-Means Clustering Tuning Parameter
#' @return The function returns a list of the following components:
#' \item{OptimalK}{The optimal number of clusters chosen from \code{nClusters}.}
#' \item{plotinfo}{A list containing the information needed in S3 method \code{plot}.}
#' @family sparse weighted K-Means functions
#' @author Wenyu Zhang
#' @references Robert, T. \emph{et al.} (2001). Estimating the number of clusters in a data set via the gap statistic. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{63(2)}, 411-423.
#' @export

ChooseK <- function(x, nClusters=2:6, nperms=20,nstart=20){
  ObjPerm <- matrix(NA,nrow = length(nClusters),ncol = nperms)
  ObjData <- numeric(length = length(nClusters))
  for(i in 1:nperms){
    permx <- matrix(NA, nrow=nrow(x), ncol=ncol(x))
    for (p in 1:ncol(x)) permx[,p] <- sample(x[,p])
    for (j in 1:length(nClusters)){
      K <- nClusters[j]
      ObjPerm[j,i] <- kmeans(permx,K,nstart = nstart)$tot.withinss
      ObjData[j] <- kmeans(x,K,nstart = nstart)$tot.withinss
    }
  }
  ObjPerm <- t(ObjPerm)
  ymean <- apply(log(ObjPerm), 2, mean) - log(ObjData)
  ysd <- apply(log(ObjPerm), 2, sd)
  BestK.pos <- which.max(ymean>=(max(ymean)-ysd[which.max(ymean)]))
  out <- list(OptimalK=nClusters[BestK.pos],plotinfo=list(nClusters=nClusters,ymean=ymean,ysd=ysd))
  class(out) <- "ChooseK"
  return(out)
}

#' @describeIn ChooseK plot the Gap statistic of each candicate number of clusters K.
#' @export
plot.ChooseK <- function(x, ...){
  nClusters <- x$plotinfo$nClusters
  ymean <- x$plotinfo$ymean
  ysd <- x$plotinfo$ysd
  plot(nClusters, ymean, ylim = range(c(ymean - ysd, ymean + 
                                          ysd)), pch = 19, xlab = "candidate number of clusters", ylab = "Gap", 
       cex = 0.3, main = paste("Optimal K =", x$OptimalK))
  lines(nClusters, ymean, col = "red", lty = 2)
  arrows(nClusters, ymean - ysd, nClusters, ymean + ysd, length = 0.05, 
         angle = 90, code = 3)
}
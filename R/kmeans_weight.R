#' Weighted K-Means Clustering with Weights on Observations
#'
#' Perform K-Means algorithm on observations with given weights.
#' @param x An \emph{n} by \emph{p} numeric data matrix, and \emph{n} is the number of observations and \emph{p} the number of features.
#' @param K The number of clusters. Omitted if \code{centers} are provided.
#' @param weight A vector of \emph{n} positive elements representing weights on observations.
#' @param centers A \emph{K} by \emph{p} matrix indicating initial (distinct) cluster centers.
#' @param nstart The number of initial random sets chosen from (distinct) rows in \code{x}. Omitted if \code{centers} is provided. Default is 20.
#' @param algorithm Character; either "\code{Hartigan-Wong}" or "\code{Forgy}". Default is "\code{Hartigan-Wong}".
#' @keywords Weighted K-Means Clustering
#' @return The function returns a list of the following components:
#' \item{centers}{the centers of the clustering result.}
#' \item{cluster}{a vector of integers (from \code{1:k}) indicating the cluster to which each observation is allocated.}
#' \item{weight}{a vector of non-zero weights in the input vector \code{weight}.}
#' \item{wcss}{normalized within-cluster sum of squares, i.e. the objective divided by \code{sum(weight)}.}
#' @family sparse weighted K-Means functions
#' @author Wenyu Zhang
#' @importFrom stats na.omit rnorm sd
#' @importFrom graphics lines arrows plot
#' @importFrom Rcpp evalCpp
#' @useDynLib SWKM
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' data("DMdata")
#' # data preprocessing
#' data <- t(DMdata$data)
#' data_rank <- apply(data, 2, rank) 
#' data_rank_center<- t(t(data_rank) - colMeans(data_rank)) 
#' data_rank_center_scale <- t(t(data_rank_center)/apply(data_rank_center, 2, sd)) 
#' data_processed <-  t(data_rank_center_scale) 
#' # tune the number of cluster K
#' # nperms and nstart are set to be small in order to save computation time
#' cK <- ChooseK(data_processed[-DMdata$noisy.label,],nClusters = 1:6,nperms = 10,nstart = 5)
#' plot(cK)
#' K <- cK$OptimalK
#' # tune weight
#'   res.tuneU <- kmeans.weight.tune(x = data_processed,K = K,
#'   noisy.lab = DMdata$noisy.label,nperms = 10,nstart = 5)
#' plot(res.tuneU)
#' # perform weighted K-means
#' res <- kmeans.weight(x = data_processed,K = K,weight = res.tuneU$bestweight)
#' # check the result
#' table(res$cluster,DMdata$true.label)
#' }

kmeans.weight <- function(x,K=NULL,weight=NULL,centers=NULL,nstart=20,algorithm="Hartigan-Wong"){
  if (is.null(K) && is.null(centers))
    stop("Must provide either K or centers.")
  if (!is.null(K) && !is.null(centers) && nrow(centers) != K )
    stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
  if (!is.null(centers) && ncol(centers)!=ncol(x))
    stop("If centers is provided, then ncol(centers) must equal ncol(x).")
  if(!is.null(weight) && length(weight)!=nrow(x))
    stop("length(weight) must equal nrow(x).")
  if(is.null(weight))    weight <- rep(1,nrow(x))
  if(any(weight<0)) stop("weight should be larger than or equal to 0.")
  #### if there exists some zero weight, the corresponding data will be omitted. ####
  if(any(weight==0)) {
    x <- x[weight!=0,]
    weight <- weight[weight!=0]
  }

  if(is.null(centers))
    centers <- x[sample(nrow(x),K),]
  else {
    nstart <- 1
    K<- nrow(centers)
  }
  if(length(weight)<=K)
    stop("number of positive weights should be larger than K and nrow(centers)!")
  if(K==1){
    centers <- weight/sum(weight)%*%x
    wcss <- sum(weight/sum(weight)%*%(sweep(x,2,centers)^2))
    cluster <- rep(1,nrow(x))
    return(list(centers=centers,cluster=cluster,weight=weight,wcss=wcss))
  }
  WCSS.min <- Inf
  for(ns in 1:nstart){
    if(ns > 1) centers <- x[sample(nrow(x),K),]
    niter <- 0
    maxiter <- 10
    #### Forgy algorithm ####
    if(algorithm=="Forgy"){
      centers.old <- centers+1
      while(niter<maxiter && sum(abs(centers-centers.old))>1e-4){
        distmat2 <- apply(centers, 1, function(y){colSums((t(x) - y)^2)})
        Cs <- apply(distmat2, 1, which.min)
        centers.old <- centers
        if(length(unique(Cs))!=K)
          stop("empty cluster: try a better set of initial centers.")
        for(k in unique(Cs))
          centers[k,] <- (weight[Cs==k]/sum(weight[Cs==k]))%*%x[Cs==k,,drop=FALSE]
        niter <- niter + 1
      }
      WCSS <- GetWCSS.weight(x,Cs,weight/sum(weight))$wcss
      if(WCSS<WCSS.min) {
        WCSS.min <- WCSS
        centers.min <- centers
        Cs.min <- Cs
      }
    }
    #### Hartigan Algorithm ####
    if(algorithm=="Hartigan-Wong"){
      out <- kmeans_hartigan(x,centers,weight/sum(weight))
      WCSS <- out$wcss
      if(WCSS<WCSS.min) {
        WCSS.min <- WCSS
        centers.min <- out$centers
        Cs.min <- out$cluster
      }
    }
  }
  return(list(centers=centers.min,cluster=Cs.min,weight=weight,wcss=WCSS.min))
}

GetWCSS.weight <- function(x, Cs, weight, ws = NULL){
  wcss.perfeature <- numeric(ncol(x))
  for (k in unique(Cs)) {
    whichers <- (Cs == k)
    if (sum(whichers) > 1) {
      centers <- (weight[whichers]/sum(weight[whichers]))%*%x[whichers,]
      wcss.perfeature <- wcss.perfeature + weight[whichers]%*%(sweep(x[whichers,],2,centers)^2)
    }
  }
  centers <- (weight/sum(weight))%*%x
  bcss.perfeature <- weight%*%(sweep(x,2,centers)^2) - wcss.perfeature
  if (!is.null(ws))
    return(list(wcss.perfeature = wcss.perfeature, wcss.ws = sum(wcss.perfeature * ws),
                bcss.perfeature = bcss.perfeature, bcss.ws = sum(bcss.perfeature * ws)))
  if (is.null(ws))
    return(list(wcss.perfeature = wcss.perfeature, wcss = sum(wcss.perfeature),
                bcss.perfeature = bcss.perfeature, bcss = sum(bcss.perfeature)))
}

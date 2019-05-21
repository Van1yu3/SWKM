#' Weighted K-Means Clustering with Weights on Observations
#'
#' Perform K-Means algorithm on observations with given weights.
#' @param x An \emph{n} by \emph{p} numeric data matrix, and \emph{n} is the number of observations and \emph{p} the number of features.
#' @param K The number of clusters. Omitted if \code{centers} are provided.
#' @param weight A vector of \emph{n} positive elements representing weights on observations.
#' @param centers A \emph{K} by \emph{p} matrix indicating initial (distinct) cluster centers.
#' @param nstart The number of initial random sets chosen from (distinct) rows in \code{x}. Omitted if \code{centers} is provided. Default is 20.
#' @param algorithm Character; either "\code{Hartigon-Wong}" or "\code{Forgy}". Default is "\code{Hartigon-Wong}".
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
#' # generate data
#' set.seed(1)
#' require(mvtnorm)
#' n <- 60  #sample size
#' p <- 1000 #dimension of features
#' q <- 50  #dimension of cluster-specific features
#' mu <- 0.8
#' MU <- c(0,-mu,mu)
#' sigma0 <- 5
#' data <- rbind(rmvnorm(n/3,rep(0,p)),rmvnorm(n/3,c(rep(-mu,q),rep(0,p-q))),
#' rmvnorm(n/3,c(rep(mu,q),rep(0,p-q))))
#' # add noise to 10 random observations
#' noisy.lab <- sample(n,10)
#' for (k in 1:3){
#' check <- (noisy.lab<n*k/3+1) & (noisy.lab>n/3*(k-1))
#' temp.lab <- noisy.lab[check]
#' num <- length(temp.lab)
#' if(any(check))
#'   data[temp.lab,] <- rmvnorm(num,c(rep(MU[k],q),rep(0,p-q)),sigma = diag(sigma0,p))
#' }
#' # apply kmeans.weight.tune to tune weight parameter U
#' res.tuneU <- kmeans.weight.tune(data,K=3,noisy.lab=noisy.lab)
#' plot(res.tuneU)
#' # perform weighted K-Means to get the clustering result
#' res <- kmeans.weight(data,K=3,weight=res.tuneU$bestweight)

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
    }
    #### Hartigan Algorithm ####
    if(algorithm=="Hartigan-Wong"){
      return(kmeans_hartigan(x,centers,weight))
    }
    WCSS <- GetWCSS.weight(x,Cs,weight/sum(weight))$wcss
    if(WCSS<WCSS.min) {
      WCSS.min <- WCSS
      centers.min <- centers
      Cs.min <- Cs
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

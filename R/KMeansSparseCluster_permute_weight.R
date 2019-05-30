#' Choose Sparsity Parameter for Sparse Weighted K-Means Clustering
#'
#' The sparsity parameter controls the L1 bound on w, the feature weights.
#'  A permutation approach is used to select the sparsity parameter.
#' @inheritParams KMeansSparseCluster.weight
#' @inheritParams kmeans.weight.tune
#' @param nvals The number of candidate tuning parameter values. Omitted if \code{wbounds} is given.
#' @inherit sparcl::KMeansSparseCluster.permute return
#' @inherit KMeansSparseCluster.weight examples
#' @keywords Sparse Weighted K-Means Clustering Tuning Parameter
#' @family sparse weighted K-Means functions 
#' @author Wenyu Zhang
#' @references Daniela M Witten and Robert Tibshirani (2010). A framework for feature selection in clustering.  \emph{Journal of the American Statistical Association}, \bold{105(490)}, 713-726.
#' @export


KMeansSparseCluster.permute.weight <-
  function(x, K=NULL, weight=NULL, nperms=20, nstart=20, wbounds=NULL,silent=TRUE, nvals=10, centers=NULL){
    if(is.null(wbounds)) wbounds <- exp(seq(log(1.2), log(sqrt(ncol(x))*.9), len=nvals))
    if(min(wbounds) <= 1) stop("Wbounds should be greater than 1, since otherwise only one weight will be nonzero.")
    if(length(wbounds)<2) stop("Wbounds should be a vector of at least two elements.")
    # was seq(1.2, sqrt(ncol(x))*.6, len=10)
    if(is.null(K) && is.null(centers)) stop("Must provide either K or centers.")
    if(!is.null(K) && !is.null(centers)){
      if(nrow(centers)!=K) stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
      if(nrow(centers)==K) K <- NULL
    }
    if(!is.null(centers) && ncol(centers)!=ncol(x)) stop("If centers is provided, then ncol(centers) must equal ncol(x).")           
    if(!is.null(weight) && length(weight)!=nrow(x))
      stop("length(weight) must equal nrow(x).")
    if(is.null(weight))
      weight <- rep(1,nrow(x))
    out <- KMeansSparseCluster.weight(x, K,weight = weight, wbounds=wbounds,nstart = nstart,
                                      silent=silent, centers=centers)
    nnonzerows <- tots <- numeric(length(out))
    for(i in 1:length(out)){
      nnonzerows[i] <- sum(out[[i]]$ws!=0)
      tots[i] <- out[[i]]$wcss$bcss.ws
    }
    permtots <- matrix(NA, nrow=length(wbounds), ncol=nperms)
    for (i in 1:nperms){
      if(!silent) cat("Permutation ", i, "of ", nperms, fill=TRUE)
      permx <- matrix(NA, nrow=nrow(x), ncol=ncol(x))
      for (j in 1:ncol(x)) permx[,j] <- sample(x[,j])
      perm.out <- KMeansSparseCluster.weight(permx, K, weight, wbounds=wbounds,nstart = nstart,
                                               silent=silent, centers=centers)
      for(j in 1:length(perm.out)) permtots[j,i] <- perm.out[[j]]$wcss$bcss.ws
    }
    gaps <- (log(tots)-apply(log(permtots),1,mean))
    out <- list(tots=tots, permtots=permtots, nnonzerows=nnonzerows, gaps=gaps, 
                sdgaps=apply(log(permtots),1,sd), wbounds=wbounds, bestw=wbounds[which.max(gaps)])
    if(!silent) cat(fill=TRUE)
    class(out) <- "KMeansSparseCluster.permute.weight"
    return(out)
  }
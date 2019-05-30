#' Sparse Weighted K-Means Clustering with Weights on Observations
#'
#' Perform sparse weighted K-Means algorithm on observations with given weights.
#' @inheritParams kmeans.weight.tune
#' @inheritParams kmeans.weight
#' @inheritParams sparcl::KMeansSparseCluster
#' @keywords Sparse Weighted K-Means Clustering
#' @return If \code{wbounds} is a numeric value, then the function returns a list
#'  with elements as follows:
#'  \item{ws}{The p-vector of feature weights.}
#'  \item{Cs}{The clustering result.}
#'  \item{wcss}{A list of the following:
#'  \code{wcss.perfeature}, \code{wcss.ws}, \code{bcss.perfeature}, \code{bcss.ws}.
#'  Among them, \code{wcss.ws}=\code{sum(wcss.perfeature*ws)},
#'  \code{bcss.ws}=\code{sum(bcss.perfeature*ws)}. And \code{bcss.ws} is the objective in
#'  sparse weighted K-Means clustering algorithm.}
#'  \item{wbound}{The L1 bound in the current list.}
#'  \item{weight}{The weights on observations.}
#'  If \code{wbounds} is a vector, then the function returns a list with lists
#'  (one per element of \code{wbounds}).
#' @family sparse weighted K-Means functions
#' @author Wenyu Zhang
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' data("NormalDisData")
#' cK <- ChooseK(NormalDisData$data[-NormalDisData$noisy.label,],nClusters = 1:6)
#' plot(cK)
#' K <- cK$OptimalK
#' res.tuneU <- kmeans.weight.tune(x = NormalDisData$data,K = K,
#' noisy.lab = NormalDisData$noisy.label,weight.seq = NULL)
#' plot(res.tuneU)
#' res.tunes <- KMeansSparseCluster.permute.weight(x = NormalDisData$data,K = K,
#' weight = res.tuneU$bestweight)
#' res <- KMeansSparseCluster.weight(x = NormalDisData$data,K = K,
#' wbounds = res.tunes$bestw,weight = res.tuneU$bestweight)
#' #check the clustering result, the number of features selected and the 50 most important features 
#' table(res[[1]]$Cs,NormalDisData$true.label)
#' sum(res[[1]]$ws!=0)
#' order(res[[1]]$ws,decreasing = TRUE)[1:50]
#' }





  

KMeansSparseCluster.weight <- function (x, K = NULL, weight=NULL, wbounds = NULL, nstart = 20,
                                        silent = TRUE, maxiter = 6, centers = NULL){
  if (is.null(K) && is.null(centers))
    stop("Must provide either K or centers.")
  if (!is.null(K) && !is.null(centers)) {
    if (nrow(centers) != K)
      stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
    if (nrow(centers) == K)
      K <- NULL
  }
  if (!is.null(centers) && ncol(centers) != ncol(x))
    stop("If centers is provided, then ncol(centers) must equal ncol(x).")
  if(!is.null(weight) && length(weight)!=nrow(x))
    stop("length(weight) must equal nrow(x).")
  if(is.null(weight))
    weight <- rep(1,nrow(x))
  if (is.null(wbounds))
    wbounds <- seq(1.1, sqrt(ncol(x)), len = 20)
  if (min(wbounds) <= 1)
    stop("wbounds should be greater than 1")
  wbounds <- c(wbounds)
  out <- list()
  if (!is.null(K))
    Cs <- kmeans.weight(x, K = K, weight = weight, nstart = nstart)$cluster
  if (is.null(K))
    Cs <- kmeans.weight(x, centers = centers,weight = weight)$cluster
  for (i in 1:length(wbounds)){
    if (length(wbounds) > 1 && !silent)
      cat(i, fill = FALSE)
    ws <- rep(1/sqrt(ncol(x)), ncol(x))
    ws.old <- rnorm(ncol(x))
    niter <- 0
    while ((sum(abs(ws - ws.old))/sum(abs(ws.old))) > 1e-04 && niter < maxiter) {
      if (!silent)
        cat(niter, fill = FALSE)
      niter <- niter + 1
      ws.old <- ws
      if (!is.null(K)) {
        if (niter > 1)
          Cs <- UpdateCs.weight(x, K, ws, Cs, weight)
      } else {
        if (niter > 1)
          Cs <- UpdateCs.weight(x, nrow(centers), ws, Cs, weight)
      }
      ws <- UpdateWs.weight(x, Cs, wbounds[i], weight)
    }
    out[[i]] <- list(ws = ws, Cs = Cs, wcss = GetWCSS.weight(x, Cs, weight, ws),
                     wbound = wbounds[i],weight = weight)
  }
  if (!silent)
    cat(fill = TRUE)
  class(out) <- "KMeansSparseCluster"
  return(out)
}

UpdateCs.weight <- function (x, K, ws, Cs, weight){
  x <- x[, ws != 0, drop=FALSE]
  z <- sweep(x, 2, sqrt(ws[ws != 0]), "*")
  nrowz <- nrow(z)
  mus <- matrix(NA,max(Cs),ncol(z))
  for (k in unique(Cs))
    mus[k,] <- (weight[Cs==k]/sum(weight[Cs==k]))%*%z[Cs==k,]
  distmat <- apply(mus, 1, function(y){colSums((t(z) - y)^2)})
  nearest <- apply(distmat, 1, which.min)
  if (length(unique(nearest)) == K) {
    km <- kmeans.weight(z, weight = weight, centers = mus)
  } else {
    km <- kmeans.weight(z, K = K, weight = weight, nstart = 20)
  }
  return(km$cluster)
}

UpdateWs.weight <-  function (x, Cs, l1bound,weight){
  wcss.perfeature <- GetWCSS.weight(x, Cs=Cs, weight=weight)$wcss.perfeature
  tss.perfeature <- GetWCSS.weight(x, rep(1, nrow(x)),weight)$wcss.perfeature
  lam <- BinarySearch(-wcss.perfeature + tss.perfeature, l1bound)
  ws.unscaled <- soft(-wcss.perfeature + tss.perfeature, lam)
  return(ws.unscaled/l2n(ws.unscaled))
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

BinarySearch <- function (argu, sumabs){
  if (l2n(argu) == 0 || sum(abs(argu/l2n(argu))) <= sumabs)
    return(0)
  lam1 <- 0
  lam2 <- max(abs(argu)) - 1e-05
  iter <- 1
  while (iter <= 15 && (lam2 - lam1) > (1e-04)){
    su <- soft(argu, (lam1 + lam2)/2)
    if (sum(abs(su/l2n(su))) < sumabs){
      lam2 <- (lam1 + lam2)/2
    }
    else {
      lam1 <- (lam1 + lam2)/2
    }
    iter <- iter + 1
  }
  return((lam1 + lam2)/2)
}

soft <- function(x, d){
  return(sign(x) * pmax(0, abs(x) - d))
}

l2n <- function(x){
  return(sqrt(sum(x^2)))
}


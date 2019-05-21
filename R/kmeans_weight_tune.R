#' Choose Tuning Parameter for Weighted K-Means Clustering
#'
#' The tuning parameter controls the weight of noisy observations.
#' A permutation approach is used to select the tuning parameter.
#' @inheritParams kmeans.weight
#' @param noisy.lab A vector indicating the row positions of noisy observations. Omitted if \code{weight.seq} is provided.
#' @param weight.seq A candidate weight matrix, each row indicating one candidate weight vector. If \code{NULL},
#' the function will assign a sequence of candidate weight \code{c(1,0.8,0.5,0.2,0.1,0.08,0.05,0.02,0.01,0.005,0.001)}
#' to noisy observations by default.
#' @param nperms Number of permutations. Default is 20.
#' @param ... unused.
#' @inherit kmeans.weight examples
#' @keywords Weighted K-Means Clustering Tuning Parameter
#' @return The function returns a list of the following components:
#' \item{gaps}{The gap statistics obtained (one for each of the candicate weights tried).
#' If \eqn{O(U)} is the objective function evaluated at the tuning parameter \eqn{U},
#' and \eqn{O*(U)} is the same quantity but for the permuted weights,
#' then Gap(\eqn{U})=mean(log\eqn{(O*(U))})-log\eqn{(O(U))}.}
#' \item{sdgaps}{The standard deviation of log\eqn{(O*(U))}}
#' \item{bestweight}{The best weight chosen by this method among all the candidate weights.}
#' @family sparse weighted K-Means functions
#' @author Wenyu Zhang
#' @export


kmeans.weight.tune <- function(x,K=NULL,noisy.lab=NULL,weight.seq=NULL,nperms=20,
                               centers=NULL,nstart=20,algorithm="Hartigan-Wong"){
  if(is.null(noisy.lab) & is.null(weight.seq))
    stop("must provide either noisy.lab or weight.seq!")
  if(!is.null(noisy.lab)) {
    if (min(noisy.lab)<1 | max(noisy.lab)>nrow(x))
    stop("noisy.lab should be a vector indicating noisy observation labels, ranging from 1 to nrow(x)!")
  }
  if(is.null(weight.seq)) {
    weight.candidate <- c(1,0.8,0.5,0.2,0.1,0.08,0.05,0.02,0.01,0.005,0.001)
    weight.seq <- matrix(1,length(weight.candidate),nrow(x))
    weight.seq[,noisy.lab] <- weight.candidate
  }
  ObjPerm <- matrix(NA,nrow = nrow(weight.seq),ncol = nperms)
  ObjData <- numeric(length = nrow(weight.seq))
  for (j in 1:nrow(weight.seq)){
    weight <- weight.seq[j,]
    ObjPerm[j,] <- replicate(nperms,expr = kmeans.weight(x,K,sample(weight),nstart = nstart,algorithm=algorithm)$wcss)
    ObjData[j] <- kmeans.weight(x,K,weight,nstart = nstart,algorithm=algorithm)$wcss
  }
  ymean <- rowMeans(log(ObjPerm)) - log(ObjData)
  ysd <- apply(log(ObjPerm), 1, sd)
  Bestw.pos <- which.max(ymean>=(max(ymean)-ysd[which.max(ymean)]))
  out <- list(gaps=ymean,sdgaps=ysd,bestweight=weight.seq[Bestw.pos,],weight.seq=weight.seq,Bestw.pos=Bestw.pos)
  class(out) <- "kmeans.weight.tune"
  return(out)
}
#' @describeIn kmeans.weight.tune plot the Gap statistic of each candicate weight vector.
#' @export
plot.kmeans.weight.tune <- function(x,...){
  plot(1:nrow(x$weight.seq), x$gaps, ylim = range(c(x$gaps - x$sdgaps, x$gaps + x$sdgaps)),
       pch = 19, xlab = "candidate weight label", ylab = "Gap",
       main = paste("best weight label =",x$Bestw.pos))
  lines(1:nrow(x$weight.seq), x$gaps, col = "red")
  arrows(1:nrow(x$weight.seq), x$gaps - x$sdgaps, 1:nrow(x$weight.seq), x$gaps + x$sdgaps,
         length = 0.05, angle = 90, code = 3)
}


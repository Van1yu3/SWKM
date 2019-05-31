## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  comment = "#>"
)

## ---- eval=FALSE, results='hide'-----------------------------------------
#  library(devtools)
#  devtools::install_github("Van1yu3/SWKM")

## ---- results='hide'-----------------------------------------------------
library(SWKM)

## ------------------------------------------------------------------------
data("NormalDisData")

## ------------------------------------------------------------------------
names(NormalDisData)

## ------------------------------------------------------------------------
set.seed(1)
cK <- ChooseK(NormalDisData$data[-NormalDisData$noisy.label,],nClusters = 1:6)
plot(cK)

## ------------------------------------------------------------------------
K <- cK$OptimalK
system.time(
  res.tuneU <- kmeans.weight.tune(x = NormalDisData$data,K = K,noisy.lab = NormalDisData$noisy.label,weight.seq = NULL)
)

## ------------------------------------------------------------------------
plot(res.tuneU)

## ------------------------------------------------------------------------
res.tunes <- KMeansSparseCluster.permute.weight(x = NormalDisData$data,K = K,weight = res.tuneU$bestweight)

## ------------------------------------------------------------------------
res <- KMeansSparseCluster.weight(x = NormalDisData$data,K = K,wbounds = res.tunes$bestw,weight = res.tuneU$bestweight)

## ------------------------------------------------------------------------
table(res[[1]]$Cs,NormalDisData$true.label)
sum(res[[1]]$ws!=0)
order(res[[1]]$ws,decreasing = TRUE)[1:50]

## ------------------------------------------------------------------------
set.seed(1)
data("DMdata")
# data preprocessing
data <- t(DMdata$data)
data_rank <- apply(data, 2, rank) 
data_rank_center<- t(t(data_rank) - colMeans(data_rank)) 
data_rank_center_scale <- t(t(data_rank_center)/apply(data_rank_center, 2, sd)) 
data_processed <-  t(data_rank_center_scale) 
# tune the number of cluster K
# nperms and nstart are set to be small in order to save computation time
cK <- ChooseK(data_processed[-DMdata$noisy.label,],nClusters = 1:6,nperms = 10,nstart = 5)
plot(cK)
K <- cK$OptimalK
# tune weight
# the process of tuning weight will take some time
system.time(
  res.tuneU <- kmeans.weight.tune(x = data_processed,K = K,noisy.lab = DMdata$noisy.label,nperms = 10,nstart = 5)
)
plot(res.tuneU)
# perform weighted K-means
res <- kmeans.weight(x = data_processed,K = K,weight = res.tuneU$bestweight)
# check the result
table(res$cluster,DMdata$true.label)


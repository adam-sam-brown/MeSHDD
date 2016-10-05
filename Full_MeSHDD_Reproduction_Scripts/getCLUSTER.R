#################################################
## getCLUSTER.R                                ##
## Use output of getDISTANCE.R and clusterboot ##
## to define robust clusters and select the    ##
## best value for k                            ##
#################################################

## Load
load('data/DIST.RData')
library(Hmisc)
library(fpc)

## Bootstrap Clusters
cluster.jaccard <- list()
for (i in 10:50) { # k between 10 and 50 inclusive
  this.boot <- clusterboot(DIST, B=100, bootmethod='boot', clustermethod=kmeansCBI, krange=i, seed=15555) # Bootstrap clusters
  cluster.jaccard[[i]] <- this.boot$bootresult # Jaccard similarities
  if (i%%10 == 0) {save(cluster.jaccard, file='data/CLUST.RData')}
}
save(cluster.jaccard, file='data/CLUST.RData')

## Best k clustering
best.boot <- clusterboot(DIST, B=10000, bootmethod='boot', clustermethod=kmeansCBI, krange=33, seed=15555)
save(best.boot, file='data/CLUSTBEST.RData')

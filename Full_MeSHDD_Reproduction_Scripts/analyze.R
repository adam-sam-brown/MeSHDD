##############################
## analyze.R                ##
## Reproduce paper analyses ##
##############################

########
# Load #
########

load('data/INDICATION.RData')
load('data/DIST.RData')
load('data/CLUSTENR.RData')
load('data/ENR.RData')
load('data/CLUST.RData')
library(ordinal)
library(Hmisc)

#############
# Figure S1 #
#############

jmeans <- unlist(lapply(cluster.jaccard[10:50],function(x) mean(as.numeric(x)))) # Get mean Jaccard Index
jse <- unlist(lapply(cluster.jaccard[10:50],function(x) mean(as.numeric(x))/sqrt(length(as.numeric(x))))) # Get Jaccard Index stanard error
errbar(x=10:50, y=jmeans, yplus=jmeans+jse, yminus=jmeans-jse, xlab='k (Number of Clusters)', ylab='Mean Jaccard Index (+/- SE)') # Plotting

###########################################################
# Drug-Drug Similarity is Predictive of Shared Indication #
###########################################################

## Preprocessing
DF <- data.frame(ind=as.vector(IND), dist=as.vector(DIST.TTD))
DF$ind <- factor(DF$ind, ordered=T)

## Ordinal Logit Model
OLM <-  clm(ind ~ dist, data=DF)


###########
# Example #
###########

## Figure 1 plotting
plot(hclust(as.dist(as.matrix(DIST)[clustlist[[20]], clustlist[[20]]])))

## Characteristics
median(unlist(lapply(clustlist, length))) # 31
sigind <- rep(NA,33)
for (i in 1:33) sigind[i] <- sum(CLUST.ENR[i,] <= 0.05)
sigindn <- list()
for (i in 1:33) sigindn[[i]] <- colnames(CLUST.ENR)[CLUST.ENR[i,] <= 0.05]
data.frame(c(1:33), unlist(lapply(clustlist, length)),sigind)

## Cluster 20 table
ENR <- ENR*nrow(ENR)*ncol(ENR)
park <- rep(NA, length(clustlist[[20]]))
schiz <- rep(NA, length(clustlist[[20]]))
for (i in 1:length(clustlist[[20]])) {
  park[i] <- min(ENR[clustlist[[20]][i], c(2345, 2923, 10318)])
  schiz[i] <- min(ENR[clustlist[[20]][i], c(506,  6050,  6092,  8691, 12276)])
}
c20df <- data.frame(drug=clustlist[[20]], schiz, park)
c20df$park <- ifelse(c20df$park >= 0.05, 0, 1)
c20df$schiz <- ifelse(c20df$schiz >= 0.05, 0, 1)


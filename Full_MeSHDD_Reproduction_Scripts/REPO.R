#################################################
## REPO.R                                      ##
## Combine indications with best clustering    ##
## to predict new indications for clusters     ##
#################################################


## Load
load('data/CLUSTBEST.RData')
load('data/INDICATION.RData')

## Parse clustering
clustlist <- list()
for (i in 1:33) {
  clustlist[[i]] <- names(which(best.boot$partition == i))
}

## Mapping from CLUST to TTD
map <- data.frame(name = names(best.boot$partition), ind = rep(NA, length(best.boot$partition)), stringsAsFactors=F)
for (i in 1:length(best.boot$partition)) {
   index <- grep(paste0('^',map$name[i],'$'), row.names(IND), ignore.case=T)
   if (length(index) > 0) map$ind[i] <- index
}
map <- map[complete.cases(map),]

## Build indication matrix
CLUST.IND <- matrix(0,33,ncol(IND.RAW),dimnames=list(c(1:33), colnames(IND.RAW))) # Initialize
for (i in 1:33) {
  this.cluster <- clustlist[[i]] # Get cluster
  this.rows <- rep(NA, length(this.cluster))
  for (j in 1:length(this.cluster)) if (this.cluster[j] %in% map$name) {this.rows[j] <- map$ind[map$name == this.cluster[j]]}
  this.rows <- this.rows[!is.na(this.rows)]
  CLUST.IND[i,] <- colSums(IND.RAW[this.rows,]) # Sum across cluster
}

## Calculate enrichment
CLUST.ENR <- matrix(NA, nrow(CLUST.IND), ncol(CLUST.IND), dimnames=list(row.names(CLUST.IND), colnames(CLUST.IND))) # Initialize
C.row <- rowSums(CLUST.IND) # Optimization
C.col <- colSums(CLUST.IND) # Optimization
C.sum <- sum(CLUST.IND) # Optimization
for (i in 1:nrow(CLUST.ENR)) { # Loop over clusters, then over indications
  rowsum <- C.row[i] # Recall
  for (j in 1:ncol(CLUST.ENR)) {
    colsum <- C.col[i] # Recall
    p.value <- phyper(CLUST.IND[i,j]-1,colsum,C.sum-colsum,rowsum,lower.tail=F) # Hypergeomtric test
    CLUST.ENR[i,j] <- p.value
  }
}

## Bonferonni correction
CLUST.ENR <- CLUST.ENR*nrow(CLUST.ENR)*ncol(CLUST.ENR)
save(CLUST.ENR, clustlist, file='data/CLUSTENR.RData')
#

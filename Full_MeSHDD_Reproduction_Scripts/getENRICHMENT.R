#################################################
## getENRICHMENT.R                             ##
## Use output of getFREQUENCY.R to calculate   ##
## enrichment of MeSH terms by chemical        ##
#################################################

## Load
load('data/FREQ.RData')

## Enrichment
ENR <- matrix(NA, nrow(FREQ), ncol(FREQ), dimnames = list(row.names(FREQ), colnames(FREQ))) # Initialize Matrix
FREQ.row <- rowSums(FREQ) # Optimization
FREQ.col <- colSums(FREQ) # Optimization
FREQ.sum <- sum(FREQ) # Optimization
for (i in 1:nrow(ENR)) { # Loop over drugs then over MeSH terms
  rowsum <- FREQ.row[i] # Recall
  for (j in 1:ncol(ENR)) {
    colsum <- FREQ.col[j] # Recall
    p.value <- phyper(FREQ[i,j]-1,colsum,FREQ.sum-colsum,rowsum,lower.tail=F) # Hypergeometric test for enrichment
    ENR[i,j] <- p.value #
  }
}

save(ENR, file='data/ENR.RData')

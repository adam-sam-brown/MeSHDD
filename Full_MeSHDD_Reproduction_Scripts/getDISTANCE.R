#################################################
## getDISTANCE.R                               ##
## Use output of getENRICHMENT.R to calculate  ##
## drug-drug distances                         ##
#################################################

## Load
load('data/ENR.RData')

## Make enrichment matrix binary
ENR <- ENR * nrow(ENR) * ncol(ENR) # Bonferroni
ENR <-  ifelse(ENR > 0.05, 0, 1) # Non-significant is 'off', Significant is 'on'

## Get distances
DIST <- dist(ENR,method='binary')

## Save
save(DIST,file='data/DIST.RData')

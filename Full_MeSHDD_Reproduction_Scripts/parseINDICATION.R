########################################
## parseINDICATION.R                  ##
## Extract drug indications from TTD  ##
########################################

## Load
load('data/DIST.RData')
DIST <- as.matrix(DIST)
library(data.table)

## Get TTD indications
TTD <- fread('raw/drug-disease_TTD2016.txt')
TTD.u <- unique(TTD$LNM) # Get unique drugs
TTD.i <- unique(unlist(strsplit(TTD$Indication,'; '))) # Get unique indications

## Compatibility
TTD.c <- TTD.u[unlist(sapply(row.names(DIST), function(x) grep(paste0('^',x,'$'), TTD.u, ignore.case=T)))] # Find compatible names

## Construct raw indication matrix
IND.RAW <- matrix(0, length(TTD.c), length(TTD.i), dimnames = list(TTD.c, TTD.i)) # Initialize raw indication matrix
for (i in 1:nrow(IND.RAW)) { # Loop over compatible TTD drugs
  IND.RAW[i,match(unlist(TTD[LNM == TTD.c[i],Indication]),TTD.i)] <- 1 # Set indications to 1
}

## Find pairwise overlap
IND <- matrix(0, length(TTD.c), length(TTD.c), dimnames = list(TTD.c, TTD.c)) # Initialize indication matrix
for (i in 1:nrow(IND)) { # Loop over drugs twice
  for (j in 1:ncol(IND)) {
    IND[i,j] <- length(intersect(which(IND.RAW[i,] == 1), which(IND.RAW[j,] == 1))) # Find overlapping terms
  }
}

## Get TTD-compatible DIST
TTD.d <- unlist(sapply(TTD.c, function(x) grep(paste0('^',x,'$'), row.names(DIST), ignore.case=T))) # Get compatible names
DIST.TTD <- DIST[TTD.d,TTD.d] # Get compatible row/cols

save(DIST.TTD, IND, IND.RAW, file='data/INDICATION.RData')

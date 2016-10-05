#################################################
## getFREQUENCY.R                              ##
## Combine output of parseMEDLINE.R to get     ##
## co-occurrence of MeSH terms and chemicals   ##
#################################################

## Load
load('data/MH2013.RData')
load('data/FDA2013.RData')

## Preprocessing for memory
uPMID <- unique(FDA2013$PMID)
MHFDA2013 <- subset(MH2013, PMID %in% uPMID, select = c('MH_term', 'PMID')) # Subset for memory

## Create Dictionaries
PMID.dict <- as.data.table(MHFDA2013$PMID)[, list(list(.I)), by = MHFDA2013$PMID] # Find indices of occurence for all PMIDs in MHFDA2013
DRUG.dict <- as.data.table(FDA2013$SubstanceName)[,list(list(.I)), by = FDA2013$SubstanceName] # Find indices of occurence for all DRUGs in FDA2013

## Get PMIDs for each drug
DRUG.PMID <- lapply(DRUG.dict$V1, function(x) unique(FDA2013$PMID[x])) # Find PMIDs for each DRUG
names(DRUG.PMID) <- DRUG.dict$FDA2013

## Get indices of MHFDA corresponding to each drug
DRUG.INDEX <- lapply(DRUG.PMID, function(x) PMID.dict[MHFDA2013 %in% x, unlist(V1)])
names(DRUG.INDEX) <- DRUG.dict$FDA2013

## Initialize Frequency Matrix
uMH <- unique(MH2013$MH) # Get all MH terms
FREQ <- matrix(0, length(DRUG.INDEX), length(uMH), dimnames = list(names(DRUG.INDEX), uMH)) # Initialize frequency matrix

## Construct Frequency Matrix
for (i in 1:length(DRUG.INDEX)) { # Loop over DRUG.INDEX
  INDICES <- DRUG.INDEX[[i]]
  # 1. Get MH terms from MHFDA2013 for INDICES
  # 2. Combine and unlist
  # 3. Tabulate
  # 4. Reorder to reflect FREQ
  COUNTS <- table(unname(unlist(MHFDA2013[INDICES, 'MH_term', with=F])))[uMH]
  FREQ[i,] <- COUNTS
}
FREQ[is.na(FREQ)] <- 0
save(FREQ, file='data/FREQ.RData')

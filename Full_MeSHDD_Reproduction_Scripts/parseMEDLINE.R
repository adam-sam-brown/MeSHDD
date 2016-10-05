#################################################
## parseMEDLINE.R                              ##
## Extract MeSH terms and chemical items from  ##
## raw MEDLINE batch files                     ##
#################################################

## Full_MH_SH_items and Chemical_items for 2013 downloaded from https://mbr.nlm.nih.gov/Download/ on Jan 18, 2016
# FDA approved drugs were downloaded from http://www.drugbank.ca/system/downloads/current/drug_links.csv.zip on Jan 18, 2016

## Libraries
library(data.table)

## Parse MeSH Headers and Save
MH2013 <- fread('raw/Full_MH_SH_items')
setnames(MH2013, c('MH_term','SH_term/null','MajorTopic','PMID','DateCreated','DateCompleted','DateRevised'))
save(MH2013, file='data/MH2013.RData')

## Parse Chemical Items and Save
CI2013 <- fread('raw/Chemical_items')
setnames(CI2013, c('SubstanceName','PMID','DateCreated','DateCompleted','DateRevised'))
save(CI2013, file='data/CI2013.RData')

## Parse FDA drugs
DB <- read.csv('raw/drug_links.csv', stringsAsFactors=F)

## Get drugs with MeSH terms
CI2013.unique <- unique(CI2013$SubstanceName) # Get unique CI SubstanceNames
CI2013.drugs <- CI2013.unique[unlist(sapply(DB$Name, function(x) grep(paste0('^',x,'$'), CI2013.unique, ignore.case=T)))] # Get proper SubstanceNames
FDA2013 <- subset(CI2013, SubstanceName %in% CI2013.drugs) # Subset for memory
save(FDA2013, file='data/FDA2013.RData')

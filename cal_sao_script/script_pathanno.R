#######################KEGG pathway sao#####################
setwd('/extDisk2/cal_sao/kallisto_results/')

library('KEGGAPI')
library('readr')

saores <- read_csv('SAO_DEG_whole_k.csv')

saoPathRaw <- getKEGGPathGenes('sao')
saoPathRaw <- sapply(saoPathRaw, function(x) {
  eachID <- sapply(strsplit(x, split = ':', fixed = TRUE), '[[', 2)
  return(eachID)
})

saoKEGG <- lapply(saoPathRaw, function(x) {
  return(x[x %in% saores$ID])
})

save(saoKEGG, file = 'saoKEGG.RData')
########################################################

#########################BioCyc pathway sao##############
library('KEGGAPI')
library('BioCycAPI')
library('foreach')
library('doMC')
library('stringr')
library('utils')

registerDoMC(8)

##~~~~~~~~~~~~~~~~~~~~~~~~change BioCycIDs format to KEGG~~~~~~~~~~~
cycIDsRaw <- getCycGenesList('GCF_000013425')
cycIDs <- str_replace(cycIDsRaw, 'CAALFMP', 'CaalfMp')
cycIDs <- str_replace(cycIDs, 'ORF', 'orf')
cycIDs[cycIDs == 'G3B3-18'] <- 'CaalfMp08'
cycIDs[cycIDs == 'G3B3-2'] <- 'orf19.1288'
cycIDs[cycIDs == 'G3B3-75'] <- 'CaalfMp11'
cycIDs[cycIDs == 'G3B3-3'] <- 'orf19.3733'
cycIDs[cycIDs == 'G3B3-1'] <- 'orf19.5211'
## noORFs <- cycIDs[!str_detect(cycIDs, 'orf')]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#########################################################

##########################KEGG cal######################
library('KEGGAPI')
library('magrittr')
library('dplyr')
library('readr')

setwd('/extDisk2/cal_sao/kallisto_results/')
degres <- read_csv('CAL_DEG_whole_k.csv')

## KEGGID
kegganno <- getProID('cal') %>%
  as_tibble %>%
  rename(., KEGGID = V1, KEGGanno = V2)

commonID <- kegganno$KEGGID %>%
  strsplit(., split = '_', fixed = TRUE) %>%
  sapply(., function(x) {ifelse(length(x) > 1, x[2], x)}) %>%
  substring(., first = 1, last = nchar(.) - 1) %>%
  as_tibble %>%
  rename(commonID = value)

kegganno %<>% bind_cols(., commonID)


## CGD with single allele
degres %<>%
  mutate(., commonID = ID %>% strsplit(., split = '_', fixed = TRUE) %>% sapply(., function(x) {paste(x[1:2], collapse = '')}))

## merge KEGG and CGD
## kegganno 6164 X 2
## degres 6208 X 18
## mergeres 6014 X 20
mergeres <- inner_join(degres, kegganno, by = 'commonID')

## KEGG path
calPathRaw <- getKEGGPathGenes('cal')

calKEGG <- lapply(calPathRaw, function(x) {
  return(mergeres$ID[mergeres$KEGGID %in% x])
})

save(calKEGG, file = 'calKEGG.RData')
########################################################


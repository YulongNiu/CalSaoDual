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

#######################BioCyc genes sao####################
## batch download BioCyc genes
setwd('/extDisk2/cal_sao/kallisto_results/')

library('BioCycAPI') ## version 0.2.1
library('ParaMisc')
library('doParallel')
library('foreach')
library('magrittr')

saocycFolder <- 'saocycgenes'

if (!dir.exists(saocycFolder)) {
  dir.create(saocycFolder)
} else {}

cycIDsRaw <- getCycGenes('GCF_000013425')

##~~~~~~~~~~~~parallel download~~~~~~~~~~~~~~~~~
registerDoParallel(cores = 8)

cutMat <- CutSeqEqu(length(cycIDsRaw), 8)

for (j in i:ncol(cutMat)) {

  print(paste0('It is running ', j, ' in a total of ', ncol(cutMat), '.'))

  cycAnno <- foreach(i = cutMat[1, j] : cutMat[2, j]) %dopar% {
    eachCycIDs <- getCycGeneInfo(cycIDsRaw[i])
    return(eachCycIDs)
  }
  names(cycAnno) <- cycIDsRaw[cutMat[1, j] : cutMat[2, j]]

  paste0('sao', cutMat[1, j], '_', cutMat[2, j], '.RData') %>%
    file.path(saocycFolder, .) %>%
    save(cycAnno, file = ., compress = 'xz')
}

stopImplicitCluster()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################


#######################BioCyc pathways sao####################
setwd('/extDisk2/cal_sao/kallisto_results/')

library('BioCycAPI') ## version 0.2.1
library('ParaMisc')
library('doParallel')
library('foreach')
library('magrittr')

saocycFolder <- 'saocycpaths'

if (!dir.exists(saocycFolder)) {
  dir.create(saocycFolder)
} else {}

saoPath <- getCycPathway('GCF_000013425')

##~~~~~~~~~~~~parallel download~~~~~~~~~~~~~~~~~
registerDoParallel(cores = 8)

cutMat <- CutSeqEqu(nrow(saoPath), 8)

for (j in 1:ncol(cutMat)) {

  print(paste0('It is running ', j, ' in a total of ', ncol(cutMat), '.'))

  cycAnno <- foreach(i = cutMat[1, j] : cutMat[2, j]) %dopar% {
    eachCycIDs <- saoPath[i, 1] %>%
      as.character %>%
      getCycGenesfPathway
    return(eachCycIDs)
  }
  names(cycAnno) <- saoPath$pathID[cutMat[1, j] : cutMat[2, j]]

  paste0('sao', cutMat[1, j], '_', cutMat[2, j], '.RData') %>%
    file.path(saocycFolder, .) %>%
    save(cycAnno, file = ., compress = 'xz')
}

stopImplicitCluster()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################################################


#########################BioCyc pathway sao##############
setwd('/extDisk2/cal_sao/kallisto_results/')

library('KEGGAPI') ## version 0.1.7.4
library('BioCycAPI') ## version 0.2.1
library('readr')
library('magrittr')
library('foreach')
library('tibble')
library('dplyr')

saores <- read_csv('SAO_DEG_whole_k.csv')

## saocyc genes
saocycFolder <- 'saocycgenes'
saoFiles <- dir(saocycFolder, full.names = TRUE)

saoConv <- foreach(i = seq_along(saoFiles), .combine = bind_rows) %do% {
  load(saoFiles[i])
  eachConv <- tibble(cycID = names(cycAnno),
         keggID = sapply(cycAnno, function(x){
           eachname <- x$name %>%
             .[grepl('SAOUHSC', x$name)]
           return(eachname)
         }))
  return(eachConv)
}
## BioCyc Sao 2923 X 2
## saores 2633 X 17

## BioCyc pathways
saocycFolder <- 'saocycpaths'
saoFiles <- dir(saocycFolder, full.names = TRUE)

saoPath <- foreach(i = seq_along(saoFiles), .combine = append) %do% {
  load(saoFiles[i])
  eachPath <- lapply(cycAnno, '[[', 1)
  return(eachPath)
}

## use kegg ID in BioCyc paths
saoBioCyc <- lapply(saoPath, function(x) {
  saoConv %>%
    filter(cycID %in% x) %>%
    .$keggID %>%
    .[. %in% saores$ID]
})
## remove no genes
saoBioCyc %<>% .[sapply(., length) > 0]

save(saoBioCyc, file = 'saoBioCyc.RData')
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

#######################BioCyc pathways cal####################
setwd('/extDisk2/cal_cal/kallisto_results/')

library('BioCycAPI') ## version 0.2.1
library('ParaMisc')
library('doParallel')
library('foreach')
library('magrittr')

calcycFolder <- 'caocycpaths'

if (!dir.exists(calcycFolder)) {
  dir.create(calcycFolder)
} else {}

calPath <- getCycPathway('CALBI')

##~~~~~~~~~~~~parallel download~~~~~~~~~~~~~~~~~
registerDoParallel(cores = 8)

cutMat <- CutSeqEqu(nrow(calPath), 8)

for (j in 73:ncol(cutMat)) {

  print(paste0('It is running ', j, ' in a total of ', ncol(cutMat), '.'))

  cycAnno <- foreach(i = cutMat[1, j] : cutMat[2, j]) %dopar% {
    eachCycIDs <- calPath[i, 1] %>%
      as.character %>%
      getCycGenesfPathway
    return(eachCycIDs)
  }
  names(cycAnno) <- calPath$pathID[cutMat[1, j] : cutMat[2, j]]

  paste0('cal', cutMat[1, j], '_', cutMat[2, j], '.RData') %>%
    file.path(calcycFolder, .) %>%
    save(cycAnno, file = ., compress = 'xz')
}

stopImplicitCluster()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################################################


#########################BioCyc pathway cal##############
setwd('/extDisk2/cal_sao/kallisto_results/')

library('KEGGAPI') ## version 0.1.7.4
library('BioCycAPI') ## version 0.2.1
library('readr')
library('magrittr')
library('foreach')
library('tibble')
library('dplyr')

calres <- read_csv('CAL_DEG_whole_k.csv')

## calcyc genes
cycIDRaws <- tibble(CycID = getCycGenes('CALBI'))
cycID <- cycIDRaws %>%
  mutate(ORF = cycIDRaws$CycID %>%
           strsplit(split = ':', fixed = TRUE) %>%
           sapply('[[', 2) %>%
           tolower)
calID <- tibble(ID = calres$ID,
                ORF = calres$Product %>%
                  strsplit(split = ')', fixed = TRUE) %>%
                  sapply('[[', 1) %>%
                  substring(first = 2))

calConv <- inner_join(cycID, calID, by = 'ORF')

## BioCyc Cal 6093 X 2
## calres 6207 X 17
## calConv 6039 X 3

## BioCyc pathways
calcycFolder <- 'calcycpaths'
calFiles <- dir(calcycFolder, full.names = TRUE)

calPath <- foreach(i = seq_along(calFiles), .combine = append) %do% {
  load(calFiles[i])
  eachPath <- lapply(cycAnno, '[[', 1)
  return(eachPath)
}

## use kegg ID in BioCyc paths
calBioCyc <- lapply(calPath, function(x) {
  calConv %>%
    filter(CycID %in% x) %>%
    .$ID %>%
    .[. %in% calres$ID]
})
## remove no genes
calBioCyc %<>% .[sapply(., length) > 0]

save(calBioCyc, file = 'calBioCyc.RData')
#########################################################


############################test merge ec#######################
setwd('/extDisk2/cal_sao/kallisto_results/')

library('Biostrings')
library('magrittr')
library('RNASeqQuant')

calnum <- readBStringSet('cal_cdna_name.fa') %>%
  length

saonum <- readBStringSet('sao_cdna_name.fa') %>%
  length

mixec <- read.table('testpseudo/pseudoalignments.ec', stringsAsFactors = FALSE) %>%
  `[`(, 2) %>%
  RNASeqQuant:::SplitEC(.)

## cal number: 6226
## sao number: 2844

##~~~~~~~~~~~~~~~~~~~~~~test mix ec in both cal and sao~~~~~~~~~~~
hasMixEc <- sapply(mixec, function(x) {
  lenx <- length(x)
  eachLogic <- sum(x < calnum) == lenx | sum((x >= calnum) & (x <= (calnum + saonum - 1))) == lenx
  return(eachLogic)
})

sum(hasMixEc) == length(mixec)

## no mix ec in cal and sao
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~split ec/efflen/ecount~~~~~~~~~~~~~~~~
mixec <- read.table('testpseudo/pseudoalignments.ec', stringsAsFactors = FALSE) %>%
  `[`(, 2)

mixeccount <- read.table('testpseudo/pseudoalignments.tsv', stringsAsFactors = FALSE) %>%
  `[`(, 2)

mixefflen <- read.table('testquant/abundance.tsv', stringsAsFactors = FALSE, header = TRUE) %>%
  `[`(, 3)

SplitSpe <- function(ec, eccount, efflen, spenum) {

  require('RNASeqQuant')
  require('magrittr')

  ## ec of the species must be completely separated
  ## two species

  numec <- RNASeqQuant:::SplitEC(ec)
  logic1st <- sapply(numec, function(x) {
    lenx <- length(x)
    eachLogic <- sum(x < spenum[1]) == lenx
    return(eachLogic)
  })

  res <- vector('list', 2)
  res[[1]]$ec <- numec[logic1st]
  res[[1]]$ec %<>% sapply(., paste, collapse = ',')
  res[[1]]$count <- eccount[logic1st]
  res[[1]]$efflen <- efflen[1:(spenum[1])]

  res[[2]]$ec <- numec[!logic1st]
  res[[2]]$ec %<>% sapply(., function(x) {
    return(paste(x - spenum[1], collapse = ','))
  })
  res[[2]]$count <- eccount[!logic1st]
  res[[2]]$efflen <- efflen[(spenum[1] + 1) : (spenum[1] + spenum[2])]

  return(res)
}

tmp1 <- SplitSpe(mixec, mixeccount, mixefflen, c(calnum, saonum))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################

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
  RNASeqQuant:::SplitEC

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
################################################################

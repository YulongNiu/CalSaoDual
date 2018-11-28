############################sao annotation##############
library('readr')
library('dplyr')
library('Biostrings')
library('stringr')
library('magrittr')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~gff~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saogff <- read_tsv('/home/Yulong/Biotools/RefData/sao/NC_007795.gff', skip = 3, col_names = c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'))

names <- saogff %>%
  select(., 'attribute') %>%
  unlist %>%
  sapply(., function(x) {
    eachNA <- x %>%
      strsplit(x, split = ';', fixed = TRUE) %>%
      unlist

    eachN <- eachNA[1] %>%
      strsplit(., split = '=', fixed = TRUE) %>%
      unlist %>%
      `[`(., 2)

    eachA <- eachNA[2] %>%
      strsplit(., split = '=', fixed = TRUE) %>%
      unlist %>%
      `[`(., 2) %>%
      substring(., first = 2, last = nchar(.) - 1)

    return(c(eachN, eachA))
  }) %>%
  t %>%
as_tibble

colnames(names) <- c('geneNames', 'Anno')

saogff %<>%
  bind_cols(., names) %>%
  mutate(., loc = paste(start, end, sep = '..'))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~add ptt/rnt~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saoptt <- read.delim('/home/Yulong/Biotools/RefData/sao/NC_007795.ptt', skip = 2, stringsAsFactor = FALSE)
saornt <- read.delim('/home/Yulong/Biotools/RefData/sao/NC_007795.rnt', skip = 2, stringsAsFactor = FALSE)
saoannot <- rbind(saoptt, saornt)

## deal with slash
slashLogic <- saoannot[, 'Gene'] == '-'
saoannot[slashLogic, 'Gene'] <- saoannot[slashLogic, 'Synonym']

saomerge <- merge(saogff, saoannot, by.x = 'loc', by.y = 'Location', sort = FALSE)

saoanno <- saomerge[, c('Synonym', 'Gene', 'V4', 'V5', 'Product'), ]
colnames(saoanno) <- c('GeneID', 'Name', 'Start', 'End', 'Product')
saoanno$Length <- abs(saoanno$Start - saoanno$End) + 1
saoanno$COG <- saomerge$COG

save(saoanno, file = '/extDisk2/cal_sao/figures_tables/saoanno.RData')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################

#################################sao k index cdna####################
library('Biostrings')
library('stringr')
library('magrittr')

load('/extDisk2/cal_sao/figures_tables/saoanno.RData')

gff <- read.table('/home/Yulong/Biotools/RefData/sao/NC_007795.gff', skip = 3, header = FALSE, sep = '\t', stringsAsFactors = FALSE) %>%
  `[`(., , 9) %>%
  strsplit(., split = ';', fixed = TRUE) %>%
  sapply(., `[`, 1) %>%
  strsplit(., split = '=', fixed = TRUE) %>%
  sapply(., `[`, 2)


## test
sum(saoanno$Name == gff) == length(gff)

smucdna <- readBStringSet('/home/Yulong/Biotools/RefData/sao/NC_007795_cdna.fa')
names(smucdna) <- saoanno$GeneID

writeXStringSet(smucdna, '/home/Yulong/Biotools/RefData/sao/NC_007795_cdna_name.fa')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############################################################

#################################cal k index cdna####################
library('Biostrings')
library('magrittr')

calcdna <- readBStringSet('/home/Yulong/Biotools/RefData/cal/C_albicans_SC5314_A22_current_default_coding.fasta.gz')

cdnaname <- names(calcdna) %>% strsplit(., split = ' ', fixed = TRUE) %>% sapply(., '[', 1)
names(calcdna) <- cdnaname

## pseudogene
## tmp1name[!(tmp1name %in% rawAnno$id)]
## [1] "C6_04200C_A" "CR_01860W_A" "C1_13530W_A" "CR_01840W_A" "CR_01850C_A"
## [6] "CR_02900W_A" "C2_06450C_A" "C5_04450C_A"

writeXStringSet(calcdna, '/extDisk2/cal_sao/kallisto_results/cal_cdna_name.fa')
#####################################################################

##########################prepare cal gff annotation#################
library('stringr')
library('utils')
library('foreach')
library('doMC')
library('magrittr')

registerDoMC(8)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~CGD raw gff~~~~~~~~~~~~~~~~~~~~~~~~
gffPath <- '/home/Yulong/Biotools/RefData/cal/C_albicans_SC5314_A22_current_features.gff'
gffAnno <- read.delim(gffPath, comment.char = '#', header = FALSE)

noteAnno <- as.character(gffAnno[, 9])
noteAnno <- str_trim(noteAnno)

ids <- sapply(strsplit(noteAnno, split = ';', fixed = TRUE), '[[', 1)
ids <- sapply(strsplit(ids, split = '=', fixed = TRUE), '[[', 2)

##"ID=CR_02900W_A;Name=CR_02900W_A;Note=(orf19.2847.1) Pseudogene; formerly an ORF Predicted by Annotation Working Group that was subsequently removed from Assembly 20;Alias=CA2870,CR_02900W,CR_02900W_B,Contig4-2393_0004,IPF17923.2,IPF17924.1,OPT2.53f,orf19.2847.1"
ExtractNote <- function(x) {
  x <- str_extract(x, 'Note=.*')
  ## check if no "Note=.*"
  if (is.na(x)) {
    return(x)
  } else {}

  x <- substring(x, 6)

  ## split with ';'
  xList <- unlist(str_split(x, ';'))
  eqIdx <- str_detect(xList, '=')
  if (sum(eqIdx) >= 1) {
    eq1stNum <- which(eqIdx)[1]
    x <- paste(xList[1:(eq1stNum-1)], collapse = ';')
  } else {}

  return(x)
}

ExtractGeneName <- function(x) {
  x <- str_extract(x, 'Gene=.*?;');

  if (!is.na(x)) {
    xLen <- nchar(x);
    x <- substr(x, 6, xLen - 1);
  } else {}

  return(x)
}

geneNames <- foreach(i = 1:length(noteAnno), .combine = c) %dopar% {
  x <- URLdecode(noteAnno[i])
  x <- ExtractGeneName(x)

  if (is.na(x)) {
    x <- ids[i]
  } else {}

  return(x)
}

noteAnno <- foreach(i = 1:length(noteAnno), .combine = c) %dopar% {
  x <- URLdecode(noteAnno[i])
  x <- ExtractNote(x)
  return(x)
}

gffAnno <- cbind(gffAnno[, -9], ids, geneNames, noteAnno)

## extract useful info
rawAnno <- gffAnno[gffAnno[, 3] == 'gene', ]
rawAnno %<>% `[`(., , c(9, 10, 1, 4, 5, 7, 11))
rawAnno$Length <- abs(rawAnno[, 4] - rawAnno[, 5]) + 1
colnames(rawAnno) <- c('GeneID', 'Names', 'Chromosome', 'Start', 'End', 'Strand', 'Product', 'Length')
write.csv(rawAnno, '/extDisk2/cal_sao/figures_tables/CGD_rawgff_Anno.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~NCBI sao gff~~~~~~~~~~~~~~~~~~~~~~~~~
gffPath <- '/home/Yulong/Biotools/RefData/sao/Sa.gff'
gffAnno <- read.delim(gffPath, comment.char = '#', header = FALSE, stringsAsFactor = FALSE)

## CDS
## ID=cds0;Parent=gene0;Dbxref=Genbank:YP_498609.1,GeneID:3919798;Name=YP_498609.1;Note=binds to the dnaA-box as an ATP-bound complex at the origin of replication during the initiation of chromosomal replication%3B can also affect transcription of multiple genes including itself.;gbkey=CDS;gene=dnaA;product=chromosomal replication initiation protein;protein_id=YP_498609.1;transl_table=11
cdsAnno <- gffAnno[gffAnno[, 3] %in% c('CDS', 'tRNA', 'rRNA'), 9]
cdsAnnoMat <- cbind(str_extract(cdsAnno, 'Parent=\\w+;?') %>% substr(., 8, nchar(.) - 1),
                    str_extract(cdsAnno, 'Note=(.*?);') %>% substr(., 6, nchar(.) - 1) %>% sapply(., URLdecode),
                    str_extract(cdsAnno, 'product=(.*?);') %>% substr(., 9, nchar(.) - 1) %>% sapply(., URLdecode))
rownames(cdsAnnoMat) <- NULL
colnames(cdsAnnoMat) <- c('ID', 'Note', 'Product')

## gene
## ID=gene0;Dbxref=GeneID:3919798;Name=dnaA;gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=SAOUHSC_00001
geneMat <- gffAnno[gffAnno[, 3] == 'gene', ]
geneAnno <- geneMat[, 9]
geneAnnoMat <- cbind(str_extract(geneAnno, 'ID=(\\w+?);') %>% substr(., 4, nchar(.) - 1),
                    str_extract(geneAnno, 'Name=(.*?);') %>% substr(., 6, nchar(.) - 1) %>% sapply(., URLdecode),
                    str_extract(geneAnno, 'locus_tag=.*') %>% substr(., 11, nchar(.)) %>% sapply(., URLdecode))
geneMat <- cbind(geneAnnoMat, geneMat[, c(4, 5, 7)])
rownames(geneMat) <- NULL

rawAnno <- merge(geneMat, cdsAnnoMat, by.x = '1', by.y = 'ID', sort = FALSE)
rawAnno %<>% `[`(, c(1, 3, 2, 4:8))
colnames(rawAnno)[1:6] <- c('ID', 'GeneID', 'Name', 'Start', 'End', 'Strand')
rawAnno$Length <- abs(rawAnno$Start - rawAnno$End) + 1
write.csv(rawAnno, '/extDisk2/cal_sao/figures_tables/sao_rawgff_Anno.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################################################################






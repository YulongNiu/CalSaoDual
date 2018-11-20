##########################prepare gff annotation#################
library('stringr')
library('utils')
library('foreach')
library('doMC')

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
rawAnno <- gffAnno[gffAnno[, 3] == 'gene', ]
write.csv(rawAnno, '/extDisk2/cal_sao/figures_tables/CGD_rawgff_Anno.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~NCBI sao gff~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################################################################

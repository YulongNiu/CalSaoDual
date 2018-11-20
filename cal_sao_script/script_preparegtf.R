###################prepare cal gff file#####################
## http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/

############################################################



###################prepare sao gff file#####################
## http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/
setwd('/home/Yulong/Biotools/RefData/sao/')

library('readr')
library('magrittr')

gffraw <- read_tsv('Sa.gff', comment = '#', col_names = c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'))

############################################################

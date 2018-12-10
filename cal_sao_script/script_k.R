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

## test k
## mixec <- read.table('Sa_CAF2_1_1_pseudo/pseudoalignments.ec', stringsAsFactors = FALSE) %>%
##   `[`(, 2)

## mixeccount <- read.table('Sa_CAF2_1_1_pseudo/pseudoalignments.tsv', stringsAsFactors = FALSE) %>%
##   `[`(, 2)

## mixefflen <- read.table('Sa_CAF2_1_1/abundance.tsv', stringsAsFactors = FALSE, header = TRUE) %>%
##   `[`(, 3)

## eachmix <- SplitSpe(mixec, mixeccount, mixefflen, c(calnum, saonum))
## eachcal <- EM(eachmix[[1]]$efflen, eachmix[[1]]$ec, eachmix[[1]]$count, spenum = calnum)
## eachsao <- EM(eachmix[[2]]$efflen, eachmix[[2]]$ec, eachmix[[2]]$count, spenum = saonum)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###############################################################

################################load kallisto####################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~load mix ec~~~~~~~~~~~~~~~~~~
library('foreach')
library('tximport')
library('rhdf5')
library('magrittr')

## cal number: 6226
## sao number: 2844
calnum <- 6226
saonum <- 2844

wd <- '/extDisk2/cal_sao/kallisto_results'
setwd(wd)

callabel <- c('CAF2_1_1', 'CAF2_1_2', 'CAF2_1_3')
files <- file.path(wd, callabel, 'abundance.h5')
names(files) <- callabel
calk <- tximport(files, type = 'kallisto', txOut = TRUE)

saolabel <- c('Sa_1', 'Sa_2', 'Sa_3')
files <- file.path(wd, saolabel, 'abundance.h5')
names(files) <- saolabel
saok <- tximport(files, type = 'kallisto', txOut = TRUE)

calsaolabel <- c('Sa_CAF2_1_1', 'Sa_CAF2_1_2', 'Sa_CAF2_1_3')
mixcountlen <- foreach (i = seq_along(calsaolabel), .combine = cbind) %do% {

  mix <- read.table(file.path(calsaolabel[i], 'abundance.tsv'), stringsAsFactors = FALSE, header = TRUE) %>%
    `[`(, c(3, 4))

  return(mix)
}
mixcountlen %<>% as.matrix

## merge cal and sao
calk$counts %<>% cbind(., mixcountlen[1 : calnum, c(2, 4, 6)])
calk$length %<>% cbind(., mixcountlen[1 : calnum, c(1, 3, 5)])

saok$counts %<>% cbind(., mixcountlen[(calnum + 1) : (calnum + saonum), c(2, 4, 6)])
saok$length %<>% cbind(., mixcountlen[(calnum + 1) : (calnum + saonum), c(1, 3, 5)])
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~DESeq2 for cal~~~~~~~~~~~~~~~~~~~~~~~
library('DESeq2')
library('dplyr')
library('readr')
library('tibble')

calanno <- read_csv('/extDisk2/cal_sao/figures_tables/CGD_rawgff_Anno.csv')

targets <- data.frame(Group = factor(c('CAL', 'CAL', 'CAL', 'CAL_SAO', 'CAL_SAO', 'CAL_SAO')), Sample = paste0('cal', 1:6))
rownames(targets) <- paste0(targets$Group, c(1:3, 1:3))
colnames(calk$counts) <- rownames(targets)
degres <- DESeqDataSetFromTximport(calk, colData = targets, design = ~Group)

## DEGs
degres <- degres[rowSums(counts(degres)) > 1, ]
degres <- DESeq(degres)

## count transformation
rld <- rlog(degres)
vst <- varianceStabilizingTransformation(degres)
resRaw <- degres %>%
  results %T>%
  summary %>%
  as.tibble

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(rld), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as.tibble %>%
  bind_cols(.,  select(resRaw, pvalue, padj, log2FoldChange)) %>%
  inner_join(., calanno, by = 'ID') %>%
  select(., ID, Gene : Length, CAL1 : log2FoldChange) %>%
  arrange(., padj)

write.csv(res, file = 'CAL_DEG_whole_k.csv', row.names = FALSE)

## padj < 0.05 & |log2FC| > 1
## padj < 0.05 & |log2FC| > log2(1.5)
resSig <- res %>%
  filter(padj < 0.05, abs(log2FoldChange) > log2(1.5), !is.na(padj))

write.csv(resSig, file = 'CAL_DEG_FC1dot5_k.csv', row.names = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heat map~~~~~~~~~~~~~~~~~~~~~~~~~~
library('pheatmap')

## use rld
heatmapCount <- resSig %>%
  top_n(nrow(.)) %>%
  select(., Gene, CAL1 : CAL_SAO3) %>%
  as.data.frame %>%
  column_to_rownames(., var = 'Gene')

annoCol <- data.frame(Group = colData(degres)[, 1])
row.names(annoCol) <- rownames(colData(degres))
annoColor <- list(Group = c(CAL = '#00C19F', CAL_SAO = '#F8766D'))
## annoRow = data.frame(GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(30, 30, 40))))
## rownames(annoRow) <- rownames(heatmapCount)
cairo_pdf('CAL_heatmap_FC1dot5.pdf')
pheatmap(heatmapCount, annotation_col = annoCol, annotation_colors = annoColor, fontsize=12, fontsize_row=5.5, annotation_legend = TRUE)
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('directlabels')
library('ggplot2')

pca <- prcomp(t(assay(rld)))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld)))
pdf('CAL_PCA.pdf')
groupCol <- c('#F8766D', '#00C19F')
ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = groupCol) +
  geom_dl(aes(label = ID, color = Group), method = 'smart.grid')
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~volcano plot~~~~~~~~~~~~~~~~~~~~~~~~~~
library('ggplot2')
library('latex2exp')

## separate up and down
sig <- res %>%
  transmute(case_when(padj < 0.05 & log2FoldChange > log2(2) & !is.na(padj) ~ 'up',
                      padj < 0.05 & log2FoldChange < -log2(2) & !is.na(padj) ~ 'down',
                      TRUE ~ 'no')) %>%
  unlist %>%
  unname

voldt <- data.frame(padj = -log10(res$padj),
                    FC = res$log2FoldChange,
                    Type = sig)
## remove padj NA and no Inf
voldt <- voldt[!(is.na(voldt$padj) | is.infinite(voldt$padj)), ]

cairo_pdf('CAL_DEG_FC1dot5_volplot.pdf')
ggplot(voldt, aes(x = FC, y = padj, colour = Type)) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values=c('forestgreen', 'grey60', 'firebrick'),
                     labels = c('down-regulate DEGs', 'unchanged genes', 'up-regulated DEGs')) +
  geom_vline(xintercept = -log2(2), linetype = 'dashed', color = 'grey70') +
  geom_vline(xintercept= log2(2), linetype = 'dashed', color = 'grey70') +
  xlim(-5, 5) +
  xlab(TeX('$\\log_{2}$(FoldChange)')) +
  ylab(TeX('$-\\log_{10}$(adjusted P-value)'))
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~DESeq2 for sao~~~~~~~~~~~~~~~~~~~~~~~
saoanno <- read_csv('/extDisk2/cal_sao/figures_tables/saoanno.csv')

targets <- data.frame(Group = factor(c('SAO', 'SAO', 'SAO', 'SAO_CAL', 'SAO_CAL', 'SAO_CAL')), Sample = paste0('sao', 1:6))
rownames(targets) <- paste0(targets$Group, c(1:3, 1:3))
colnames(saok$counts) <- rownames(targets)
degres <- DESeqDataSetFromTximport(saok, colData = targets, design = ~Group)


## DEGs
degres <- degres[rowSums(counts(degres)) > 1, ]
degres <- DESeq(degres)

## count transformation
rld <- rlog(degres)
vst <- varianceStabilizingTransformation(degres)
resRaw <- degres %>%
  results %T>%
  summary %>%
  as.tibble

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(rld), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as.tibble %>%
  bind_cols(.,  select(resRaw, pvalue, padj, log2FoldChange)) %>%
  inner_join(., saoanno, by = 'ID') %>%
  select(., ID, Gene : Length, SAO1 : log2FoldChange) %>%
  arrange(., padj)

write.csv(res, file = 'SAO_DEG_whole_k.csv', row.names = FALSE)

## padj < 0.05 & |log2FC| > 1
## padj < 0.05 & |log2FC| > log2(1.5)
resSig <- res %>%
  filter(padj < 0.05, abs(log2FoldChange) > log2(1.5), !is.na(padj))

write.csv(resSig, file = 'SAO_DEG_FC1dot5_k.csv', row.names = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heat map~~~~~~~~~~~~~~~~~~~~~~~~~~
library('pheatmap')

## use rld
heatmapCount <- resSig %>%
  top_n(50) %>%
  select(., Gene, SAO1 : SAO_CAL3) %>%
  as.data.frame %>%
  column_to_rownames(., var = 'Gene')

annoCol <- data.frame(Group = colData(degres)[, 1])
row.names(annoCol) <- rownames(colData(degres))
annoColor <- list(Group = c(SAO = '#00C19F', SAO_CAL = '#F8766D'))
## annoRow = data.frame(GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(30, 30, 40))))
## rownames(annoRow) <- rownames(heatmapCount)
cairo_pdf('SAO_heatmap_FC1dot5.pdf')
pheatmap(heatmapCount, annotation_col = annoCol, annotation_colors = annoColor, fontsize=12, fontsize_row=5.5, annotation_legend = TRUE)
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('directlabels')
library('ggplot2')

pca <- prcomp(t(assay(rld)))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld)))
pdf('SAO_PCA.pdf')
groupCol <- c('#F8766D', '#00C19F')
ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = groupCol) +
  geom_dl(aes(label = ID, color = Group), method = 'smart.grid')
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~volcano plot~~~~~~~~~~~~~~~~~~~~~~~~~~
library('ggplot2')
library('latex2exp')

## separate up and down
sig <- res %>%
  transmute(case_when(padj < 0.05 & log2FoldChange > log2(2) & !is.na(padj) ~ 'up',
                      padj < 0.05 & log2FoldChange < -log2(2) & !is.na(padj) ~ 'down',
                      TRUE ~ 'no')) %>%
  unlist %>%
  unname

voldt <- data.frame(padj = -log10(res$padj),
                    FC = res$log2FoldChange,
                    Type = sig)
## remove padj NA and no Inf
voldt <- voldt[!(is.na(voldt$padj) | is.infinite(voldt$padj)), ]

cairo_pdf('SAO_DEG_FC1dot5_volplot.pdf')
ggplot(voldt, aes(x = FC, y = padj, colour = Type)) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values=c('forestgreen', 'grey60', 'firebrick'),
                     labels = c('down-regulate DEGs', 'unchanged genes', 'up-regulated DEGs')) +
  geom_vline(xintercept = -log2(2), linetype = 'dashed', color = 'grey70') +
  geom_vline(xintercept= log2(2), linetype = 'dashed', color = 'grey70') +
  xlim(-8, 8) +
  xlab(TeX('$\\log_{2}$(FoldChange)')) +
  ylab(TeX('$-\\log_{10}$(adjusted P-value)'))
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################

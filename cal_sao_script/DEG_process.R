###########################CAL####################
setwd('/extDisk2/cal_sao/figures_tables')

library('DESeq2')
library('magrittr')

calanno <- read.csv('CGD_rawgff_Anno.csv', row.names = 1, stringsAsFactors = FALSE)
calcount <- read.csv('cal_counts.txt', comment.char = '#', sep = '\t', stringsAsFactors = FALSE)
htCountSelect <- calcount[, c(10:12, 16:18)]
rownames(htCountSelect) <- calcount[, 1]
colnames(htCountSelect) <- c('CAL_1', 'CAL_2', 'CAL_3', 'CAL_SAO_1', 'CAL_SAO_2', 'CAL_SAO_3')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DESeq2 analysis~~~~~~~~~~~~~~~~~~~~
targets <- data.frame(Group = factor(c('CAL', 'CAL', 'CAL', 'CAL_SAO', 'CAL_SAO', 'CAL_SAO')), Sample = paste0('cal', 1:6))
rownames(targets) <- paste0(targets$Group, c(1:3, 1:3))
colnames(htCountSelect) <- rownames(targets)
glioPR <- DESeqDataSetFromMatrix(countData = htCountSelect, colData = targets, design = ~Group)

glioPR <- glioPR[rowSums(counts(glioPR)) > 1, ]
glioPR <- DESeq(glioPR)
## count transformation
rld <- rlog(glioPR)
vst <- varianceStabilizingTransformation(glioPR)
resRaw <- results(glioPR)
summary(resRaw)
res <- cbind(as.matrix(mcols(glioPR)[, 1:10]), assay(rld))
anno <- calanno[match(rownames(res), calanno[, 1]), ]
res <- cbind(anno, res[, 11:16], data.frame(resRaw[, c(5, 6, 2)]))
res <- res[order(res[, 'padj']), ]
write.csv(res, file = 'CAL_DEG_whole.csv', row.names = FALSE)

## padj < 0.01 & |log2FC| > 1
sigLogic <- res$padj < 0.05 & abs(res$log2FoldChange) > log2(2)
sigLogic[is.na(sigLogic)] <- FALSE
resSig <- res[sigLogic, ]
write.csv(resSig, file = 'CAL_DEG_FC2.csv', row.names = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heat map~~~~~~~~~~~~~~~~~~~~~~~~~~
library('pheatmap')

## use rld
topNum <- nrow(resSig)
heatmapCount <- resSig[1:topNum, 9:14]
heatmapCount <- apply(heatmapCount, 1:2, as.numeric)
rownames(heatmapCount) <- resSig[1:topNum, 'Names']

annoCol <- data.frame(Group = colData(glioPR)[, 1])
row.names(annoCol) <- rownames(colData(glioPR))
annoColor <- list(Group = c(CAL = '#00C19F', CAL_SAO = '#F8766D'))
## annoRow = data.frame(GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(30, 30, 40))))
## rownames(annoRow) <- rownames(heatmapCount)
cairo_pdf('CAL_heatmap.pdf')
pheatmap(heatmapCount, annotation_col = annoCol, annotation_colors = annoColor, fontsize=12, fontsize_row=4.5, annotation_legend = TRUE)
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
sig <- rep('no', nrow(res))
sig[res$log2FoldChange > 0 & sigLogic]  <- 'up'
sig[res$log2FoldChange < 0 & sigLogic]  <- 'down'

voldt <- data.frame(padj = -log10(res$padj),
                    FC = res$log2FoldChange,
                    Type = sig)
## remove padj NA and no Inf
voldt <- voldt[!(is.na(voldt$padj) | is.infinite(voldt$padj)), ]

cairo_pdf('CAL_DEG_FC2_volplot.pdf')
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
###########################################################



###########################SAO####################
setwd('/extDisk2/cal_sao/figures_tables')

library('DESeq2')
library('magrittr')

saoanno <- read.csv('sao_rawgff_Anno.csv', row.names = 1, stringsAsFactors = FALSE)
saocount <- read.csv('sao_counts.txt', comment.char = '#', sep = '\t', stringsAsFactors = FALSE)
htCountSelect <- saocount[, c(10:15)]
rownames(htCountSelect) <- saocount[, 1]
colnames(htCountSelect) <- c('SAO_1', 'SAO_2', 'SAO_3', 'SAO_CAL_1', 'SAO_CAL_2', 'SAO_CAL_3')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DESeq2 analysis~~~~~~~~~~~~~~~~~~~~
targets <- data.frame(Group = factor(c('SAO', 'SAO', 'SAO', 'SAO_CAL', 'SAO_CAL', 'SAO_CAL')), Sample = paste0('sao', 1:6))
rownames(targets) <- paste0(targets$Group, c(1:3, 1:3))
colnames(htCountSelect) <- rownames(targets)
glioPR <- DESeqDataSetFromMatrix(countData = htCountSelect, colData = targets, design = ~Group)

glioPR <- glioPR[rowSums(counts(glioPR)) > 1, ]
glioPR <- DESeq(glioPR)
## count transformation
rld <- rlog(glioPR)
vst <- varianceStabilizingTransformation(glioPR)
resRaw <- results(glioPR)
summary(resRaw)
res <- cbind(as.matrix(mcols(glioPR)[, 1:10]), assay(rld))
anno <- saoanno[match(rownames(res), saoanno[, 1]), ]
res <- cbind(anno, res[, 11:16], data.frame(resRaw[, c(5, 6, 2)]))
res <- res[order(res[, 'padj']), ]
write.csv(res, file = 'SAO_DEG_whole.csv', row.names = FALSE)

## padj < 0.01 & |log2FC| > 1
sigLogic <- res$padj < 0.05 & abs(res$log2FoldChange) > log2(2)
sigLogic[is.na(sigLogic)] <- FALSE
resSig <- res[sigLogic, ]
write.csv(resSig, file = 'SAO_DEG_FC2.csv', row.names = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heat map~~~~~~~~~~~~~~~~~~~~~~~~~~
library('pheatmap')

## use rld
topNum <- nrow(resSig)
heatmapCount <- resSig[1:topNum, 10:15]
heatmapCount <- apply(heatmapCount, 1:2, as.numeric)
rownames(heatmapCount) <- resSig[1:topNum, 'Name']

annoCol <- data.frame(Group = colData(glioPR)[, 1])
row.names(annoCol) <- rownames(colData(glioPR))
annoColor <- list(Group = c(SAO = '#00C19F', SAO_CAL = '#F8766D'))
## annoRow = data.frame(GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(30, 30, 40))))
## rownames(annoRow) <- rownames(heatmapCount)
cairo_pdf('SAO_heatmap.pdf')
pheatmap(heatmapCount, annotation_col = annoCol, annotation_colors = annoColor, fontsize=12, fontsize_row=4.5, annotation_legend = TRUE)
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~volcano plot~~~~~~~~~~~~~~~~~~~~~~~~~~
library('ggplot2')
library('latex2exp')

## separate up and down
sig <- rep('no', nrow(res))
sig[res$log2FoldChange > 0 & sigLogic]  <- 'up'
sig[res$log2FoldChange < 0 & sigLogic]  <- 'down'

voldt <- data.frame(padj = -log10(res$padj),
                    FC = res$log2FoldChange,
                    Type = sig)
## remove padj NA and no Inf
voldt <- voldt[!(is.na(voldt$padj) | is.infinite(voldt$padj)), ]

cairo_pdf('SAO_DEG_FC2_volplot.pdf')
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
###########################################################


###########################k res####################
library('DESeq2')
library('magrittr')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DESeq2 analysis~~~~~~~~~~~~~~~~~~~~
targets <- data.frame(Group = factor(c('CAL', 'CAL', 'CAL', 'CAL_SAO', 'CAL_SAO', 'CAL_SAO')), Sample = paste0('cal', 1:6))
rownames(targets) <- paste0(targets$Group, c(1:3, 1:3))
colnames(htCountSelect$counts) <- rownames(targets)
glioPR <- DESeqDataSetFromTximport(htCountSelect, colData = targets, design = ~Group)

glioPR <- glioPR[rowSums(counts(glioPR)) > 1, ]
glioPR <- DESeq(glioPR)
## count transformation
rld <- rlog(glioPR)
vst <- varianceStabilizingTransformation(glioPR)
resRaw <- results(glioPR)
resRaw[, 2] <- -resRaw[, 2]
summary(resRaw)
res <- cbind(as.matrix(mcols(glioPR)[, 1:10]), assay(rld))
anno <- deg[match(rownames(res), deg[, 1]), 1:8]
res <- cbind(anno, res[, 11:16], data.frame(resRaw[, c(5, 6, 2)]))
res <- res[order(res[, 'padj']), ]
Write.Csv(res, file = '/extDisk1/RESEARCH/smuSeqSongYing/kallisto_results/degseq24h_whole.csv', row.names = FALSE)

## padj < 0.01 & |log2FC| > 1
sigLogic <- res$padj < 0.01 & abs(res$log2FoldChange) > log2(1.5)
sigLogic[is.na(sigLogic)] <- FALSE
resSig <- res[sigLogic, ]
write.csv(resSig, file = '/extDisk1/RESEARCH/smuSeqSongYing/kallisto_results/degseq24h_DEG.csv', row.names = FALSE)
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
annoColor <- list(Group = c(WT = '#00C19F', deltasrtA = '#F8766D'))
## annoRow = data.frame(GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(30, 30, 40))))
## rownames(annoRow) <- rownames(heatmapCount)
cairo_pdf('/extDisk1/RESEARCH/smuSeqSongYing/kallisto_results/degseq24h_top50_heatmap.pdf')
pheatmap(heatmapCount, annotation_col = annoCol, annotation_colors = annoColor, fontsize=12, fontsize_row=7, annotation_legend = TRUE)
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
pdf('/extDisk1/RESEARCH/smuSeqSongYing/kallisto_results/degseq24h_PCA.pdf')
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

cairo_pdf('/extDisk1/RESEARCH/smuSeqSongYing/kallisto_results/degseq24h_DEG_FC1dot5_volplot.pdf')
ggplot(voldt, aes(x = FC, y = padj, colour = Type)) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values=c('forestgreen', 'grey60', 'firebrick'),
                     labels = c('down-regulate DEGs', 'unchanged genes', 'up-regulated DEGs')) +
  geom_vline(xintercept = -log2(1.5), linetype = 'dashed', color = 'grey70') +
  geom_vline(xintercept= log2(1.5), linetype = 'dashed', color = 'grey70') +
  xlim(-5, 5) +
  xlab(TeX('$\\log_{2}$(FoldChange)')) +
  ylab(TeX('$-\\log_{10}$(adjusted P-value)'))
dev.off()

cairo_pdf('/extDisk1/RESEARCH/smuSeqSongYing/kallisto_results/degseq24h_DEG_FC2_volplot.pdf')
ggplot(voldt, aes(x = FC, y = padj, colour = Type)) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values=c('forestgreen', 'grey60', 'firebrick'),
                     labels = c('down-regulate DEGs', 'unchanged genes', 'up-regulated DEGs')) +
  geom_vline(xintercept = -1, linetype = 'dashed', color = 'grey70') +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey70') +
  xlim(-5, 5) +
  xlab(TeX('$\\log_{2}$(FoldChange)')) +
  ylab(TeX('$-\\log_{10}$(adjusted P-value)'))
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~edgeR~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cts <- htCountSelect$counts
normMat <- htCountSelect$length
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y <- scaleOffset(y, t(t(log(normMat)) + o))

## filtering
keep <- filterByExpr(y)
y <- y[keep, ]

group <- factor(rep(1:2, each = 3))
design <- model.matrix(~group)
y <- estimateDisp(y, design)

## quasi-likelihood F-tests
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit,coef=2)
res <- topTags(qlf, n = 1949)
res <- res$table
sigLogic <- res$FDR < 0.01 & abs(res$logFC) > 1
sigLogic[is.na(sigLogic)] <- FALSE
resSig <- res[sigLogic, ]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################################################

############################compare DEG####################
library('magrittr')
library('venn')

setwd('/extDisk1/RESEARCH/smuSeqSongYing/kallisto_results')

##~~~~~~~~~~~~~~~~~~~~~~~~~~with anno~~~~~~~~~~~~~~~~~~~~~~~~~~~~
anno24h <- read.csv('/extDisk1/RESEARCH/smuSeqSongYing/targetDGE/24h_DEG_anno.txt', stringsAsFactor = FALSE, header = FALSE, row.names = NULL, sep = '\t')

anno4h <- read.csv('/extDisk1/RESEARCH/smuSeqSongYing/targetDGE/4h_DEG_anno.txt', stringsAsFactor = FALSE, header = FALSE, row.names = NULL, sep = '\t')

degseq24h <- read.csv('degseq24h_whole.csv', stringsAsFactor = FALSE)

degseq4h <- read.csv('degseq4h_whole.csv', stringsAsFactor = FALSE)

test24h <- merge(degseq24h, anno24h, by.y = 'V1', by.x = 'GeneID')
write.csv(test24h, 'test24h.csv')

test4h <- merge(degseq4h, anno4h, by.y = 'V1', by.x = 'GeneID')
write.csv(test4h, 'test4h.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~compare DEG~~~~~~~~~~~~~~~~~~~~~~~
deg4h <- read.csv('degseq4h_DEG_FC2.csv', stringsAsFactor = FALSE)
deg24h <- read.csv('degseq24h_DEG_FC2.csv', stringsAsFactor = FALSE)

colnames(deg4h)[9:17] %<>% paste(., '4h', sep = '_')
colnames(deg24h)[9:17] %<>% paste(., '24h', sep = '_')

deg24h %<>% `[`(., , c(1, 9:17))
commonDEG <- merge(deg4h, deg24h, by.x = 'GeneID', by.y = 'GeneID')

write.csv(commonDEG, 'degseq_4h24h_DEG_FC1dot5.csv')

## venn plot
cairo_pdf('venn_FC2.pdf')
venn(list(deg4h$GeneID, deg24h$GeneID),
     snames = c('Exponential phase DEGs', 'Stationary phase DEGs'),
     ilab = TRUE,
     zcolor = 'style',
     size = 25,
     cexil = 1.2,
     cexsn = 1.5)
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################################################################



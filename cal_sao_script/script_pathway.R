##########################GO analysis###########################
setwd('/extDisk2/cal_sao/kallisto_results')

library('goseq')
library('GO.db')
library('foreach')
library('doMC')
library('KEGGAPI')
library('BioCycAPI')
library('magrittr')
library('dplyr')
library('tibble')
library('readr')

registerDoMC(8)

load('saoGO.RData')
load('saoKEGG.RData')
load('saoBioCyc.RData')
res <- read_csv('SAO_DEG_whole_k.csv')

## remove 0 terms
saoGO %<>% `[`(sapply(saoGO, length) > 0)
saoKEGG %<>% `[`(sapply(saoKEGG, length) > 0)
saoBioCyc %<>% `[`(sapply(saoBioCyc, length) > 0)

## padj < 0.05 & |log2FC| > 1
## padj < 0.05 & |log2FC| > log2(1.5)
degVec <- res %>%
  transmute(padj < 0.05 & abs(log2FoldChange) > log2(2) & !is.na(padj)) %>%
  unlist %>%
  as.integer
names(degVec) <- res$ID

pwf <- nullp(degVec, bias.data = res$Length)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GOMat <- foreach(i = 1:length(saoGO), .combine = rbind) %dopar% {
  eachMat <- cbind(saoGO[[i]], names(saoGO)[i])
  return(eachMat)
} %>% as.data.frame

GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
  as.tibble %>%
  filter(!is.na(ontology))

## add ablog2FC
abLogFC <- sapply(GOTestWithCat$category, function(x) {
  eachres <- res %>%
    filter(., ID %in% saoGO[[x]]) %>%
    transmute(., log2FoldChange) %>%
    unlist %>%
    abs %>%
    mean(., na.rm = TRUE)

  return(eachres)
})

abLogFC %<>%
  as.data.frame %>%
  rownames_to_column(., var = 'category') %>%
  rename(., abLogFC = `.`) %>%
  as_tibble

GOTestWithCat %<>% inner_join(., abLogFC, by = 'category')

## deal with NA and select BP MF and CC
termCat <- c('BP', 'MF', 'CC')
for (i in termCat) {
  write.csv(GOTestWithCat %>% filter(ontology == i),
            paste0('SAO_FC2_', i, '_withcat.csv'))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~KEGG~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## deal path
pathAnno <- getKEGGPathAnno('sao') %>%
  as_tibble %>%
  mutate(Annotation = Annotation %>% substr(., 1, nchar(.) - 47))

KEGGMat <- foreach(i = seq_along(saoKEGG), .combine = rbind) %dopar% {
  eachMat <- cbind(saoKEGG[[i]], names(saoKEGG)[i])
  return(eachMat)
} %>% as.data.frame

KEGGTestWithCat <- goseq(pwf, gene2cat = KEGGMat, use_genes_without_cat = FALSE) %>%
  as_tibble %>%
  inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
  mutate(ontology = 'KEGG')

abLogFC <- sapply(KEGGTestWithCat$category, function(x) {
  eachres <- res %>%
    filter(., ID %in% saoKEGG[[x]]) %>%
    transmute(., log2FoldChange) %>%
    unlist %>%
    abs %>%
    mean(., na.rm = TRUE)

  return(eachres)
})

abLogFC %<>%
  as.data.frame %>%
  rownames_to_column(., var = 'category') %>%
  rename(., abLogFC = `.`) %>%
  as_tibble

KEGGTestWithCat %<>% inner_join(., abLogFC, by = 'category')

write.csv(KEGGTestWithCat, file = 'SAO_FC1dot5_KEGG_withcat.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~BioCyc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## deal path
pathAnno <- getCycPathway('GCF_000013425') %>%
  rename(Annotation = pathAnno)

BioCycMat <- foreach(i = seq_along(saoBioCyc), .combine = rbind) %dopar% {
  eachMat <- cbind(saoBioCyc[[i]], names(saoBioCyc)[i])
  return(eachMat)
} %>% as.data.frame

BioCycTestWithCat <- goseq(pwf, gene2cat = BioCycMat, use_genes_without_cat = FALSE) %>%
  as_tibble %>%
  inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
  mutate(ontology = 'BioCyc')

abLogFC <- sapply(BioCycTestWithCat$category, function(x) {
  eachres <- res %>%
    filter(., ID %in% saoBioCyc[[x]]) %>%
    transmute(., log2FoldChange) %>%
    unlist %>%
    abs %>%
    mean(., na.rm = TRUE)

  return(eachres)
})

abLogFC %<>%
  as.data.frame %>%
  rownames_to_column(., var = 'category') %>%
  rename(., abLogFC = `.`) %>%
  as_tibble

BioCycTestWithCat %<>% inner_join(., abLogFC, by = 'category')

write.csv(BioCycTestWithCat, file = 'SAO_FC2_BioCyc_withcat.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################




########################geneset analysis cal#######################
setwd('/extDisk2/cal_sao/kallisto_results')

library('goseq')
library('GO.db')
library('foreach')
library('doMC')
library('KEGGAPI')
library('BioCycAPI')
library('magrittr')
library('dplyr')
library('tibble')
library('readr')

registerDoMC(8)

load('calGO.RData')
load('calKEGG.RData')
load('calBioCyc.RData')
res <- read_csv('CAL_DEG_whole_k.csv')

## remove 0 terms
calGO %<>% `[`(sapply(calGO, length) > 0)
calKEGG %<>% `[`(sapply(calKEGG, length) > 0)
calBioCyc %<>% `[`(sapply(calBioCyc, length) > 0)

## padj < 0.05 & |log2FC| > log2(1.5)
degVec <- res %>%
  transmute(padj < 0.05 & abs(log2FoldChange) > log2(2) & !is.na(padj)) %>%
  unlist %>%
  as.integer
names(degVec) <- res$ID

pwf <- nullp(degVec, bias.data = res$Length)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GOMat <- foreach(i = 1:length(calGO), .combine = rbind) %dopar% {
  eachMat <- cbind(calGO[[i]], names(calGO)[i])
  return(eachMat)
} %>% as.data.frame

GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
  as.tibble %>%
  filter(!is.na(ontology))

## add ablog2FC
abLogFC <- sapply(GOTestWithCat$category, function(x) {
  eachres <- res %>%
    filter(., ID %in% calGO[[x]]) %>%
    transmute(., log2FoldChange) %>%
    unlist %>%
    abs %>%
    mean(., na.rm = TRUE)

  return(eachres)
})

abLogFC %<>%
  as.data.frame %>%
  rownames_to_column(., var = 'category') %>%
  rename(., abLogFC = `.`) %>%
  as_tibble

GOTestWithCat %<>% inner_join(., abLogFC, by = 'category')

## deal with NA and select BP MF and CC
termCat <- c('BP', 'MF', 'CC')
for (i in termCat) {
  write.csv(GOTestWithCat %>% filter(ontology == i),
            paste0('CAL_FC2_', i, '_withcat.csv'))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~KEGG~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## deal path
pathAnno <- getKEGGPathAnno('cal') %>%
  as_tibble %>%
  mutate(Annotation = Annotation %>% substr(., 1, nchar(.) - 19))

KEGGMat <- foreach(i = seq_along(calKEGG), .combine = rbind) %dopar% {
  eachMat <- cbind(calKEGG[[i]], names(calKEGG)[i])
  return(eachMat)
} %>% as.data.frame

KEGGTestWithCat <- goseq(pwf, gene2cat = KEGGMat, use_genes_without_cat = FALSE) %>%
  as_tibble %>%
  inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
  mutate(ontology = 'KEGG')

abLogFC <- sapply(KEGGTestWithCat$category, function(x) {
  eachres <- res %>%
    filter(., ID %in% calKEGG[[x]]) %>%
    transmute(., log2FoldChange) %>%
    unlist %>%
    abs %>%
    mean(., na.rm = TRUE)

  return(eachres)
})

abLogFC %<>%
  as.data.frame %>%
  rownames_to_column(., var = 'category') %>%
  rename(., abLogFC = `.`) %>%
  as_tibble

KEGGTestWithCat %<>% inner_join(., abLogFC, by = 'category')

write.csv(KEGGTestWithCat, file = 'CAL_FC2_KEGG_withcat.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~BioCyc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## deal path
pathAnno <- getCycPathway('CALBI') %>%
  rename(Annotation = pathAnno)

BioCycMat <- foreach(i = seq_along(calBioCyc), .combine = rbind) %dopar% {
  eachMat <- cbind(calBioCyc[[i]], names(calBioCyc)[i])
  return(eachMat)
} %>% as.data.frame

BioCycTestWithCat <- goseq(pwf, gene2cat = BioCycMat, use_genes_without_cat = FALSE) %>%
  as_tibble %>%
  inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
  mutate(ontology = 'BioCyc')

abLogFC <- sapply(BioCycTestWithCat$category, function(x) {
  eachres <- res %>%
    filter(., ID %in% calBioCyc[[x]]) %>%
    transmute(., log2FoldChange) %>%
    unlist %>%
    abs %>%
    mean(., na.rm = TRUE)

  return(eachres)
})

abLogFC %<>%
  as.data.frame %>%
  rownames_to_column(., var = 'category') %>%
  rename(., abLogFC = `.`) %>%
  as_tibble

BioCycTestWithCat %<>% inner_join(., abLogFC, by = 'category')

write.csv(BioCycTestWithCat, file = 'CAL_FC2_BioCyc_withcat.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################

##################plot pathway###########################
setwd('/extDisk2/cal_sao/kallisto_results/')

library('ggplot2')
library('RColorBrewer')
library('latex2exp')
library('magrittr')
library('stringr')

fname <- dir(pattern = 'CAL_FC2_.*_withcat.csv')
fbasename <- fname %>%
  strsplit(., split = '.', fixed = TRUE) %>%
  sapply(., `[`, 1)


for (i in seq_along(fname)) {
##for (i in 1:2) {
  pMat <- read.csv(fname[i], row.names = 1, stringsAsFactor = FALSE)

  ## pvalue < 0.05
  plotpMat <- pMat[, c('term', 'numDEInCat', 'over_represented_pvalue', 'abLogFC')]
  ## plotpMat$Annotation %<>% str_replace('<i>', '') %>% str_replace('</i>', '') %<>% str_replace('<sub>', '') %>% str_replace('</sub>', '')

  plotpMat[, 3] <- -log10(plotpMat[, 3])
  colnames(plotpMat) <- c('Name', 'Size', 'logpvalue', 'ablogFC')

  colorPal <- colorRampPalette(rev(c('red', 'yellow', 'cyan', 'blue')), bias=1)(10)

  ggplot(plotpMat[plotpMat[, 'logpvalue'] >= -log10(0.05), ], aes(x = ablogFC, y = Name)) +
  ## ggplot(plotpMat[1:10, ], aes(x = ablogFC, y = Name)) +
    geom_point(aes(size = Size, colour = logpvalue)) +
    scale_size_continuous(name = 'Number of significant genes', range = c(3,8)) +
    scale_colour_gradientn(name = '-log10(P-value)', limits=c(0, max(plotpMat[, 3])), colours = colorPal) +
    ylab('') +
    xlab(TeX('Average |$\\log_{2}$FC|')) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size = 13, face = 'bold'))
  ggsave(paste0(fbasename[i], '.pdf'), width = 13)
}
#########################################################

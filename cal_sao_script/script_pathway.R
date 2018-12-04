##########################GO analysis###########################
setwd('/extDisk2/cal_sao/kallisto_results')

library('goseq')
library('GO.db')
library('foreach')
library('doMC')
library('KEGGAPI')
library('magrittr')
library('dplyr')
library('tibble')
library('readr')

registerDoMC(8)

load('saoGO.RData')
load('saoKEGG.RData')
res <- read_csv('SAO_DEG_whole_k.csv')

## remove 0 terms
saoGO %<>% `[`(sapply(saoGO, length) > 0)
saoKEGG %<>% `[`(sapply(saoKEGG, length) > 0)

## padj < 0.05 & |log2FC| > 1
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
}
GOMat <- as.data.frame(GOMat)
GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE)
GOTestWithCat <- GOTestWithCat[!is.na(GOTestWithCat$ontology), ]

## add ablog2FC
goSub <- saoGO[match(GOTestWithCat[, 1], names(saoGO))]
abLogFC <- sapply(goSub, function(x) {
  eachFC <- res[match(x, res$GeneID), 'log2FoldChange']
  return(mean(abs(eachFC), na.rm = TRUE))
})
GOTestWithCat$abLogFC <- abLogFC

## deal with NA and select BP MF and CC
termCat <- c('BP', 'MF', 'CC')
for (i in termCat) {
  write.csv(GOTestWithCat[GOTestWithCat$ontology == i, ],
            paste0('degseq4h_FC2_', i, '_withcat.csv'))
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

write.csv(KEGGTestWithCat, file = 'SAO_FC2_KEGG_withcat.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################




########################geneset analysis cal#######################
setwd('/extDisk2/cal_sao/kallisto_results')

library('goseq')
library('GO.db')
library('foreach')
library('doMC')
library('KEGGAPI')
library('magrittr')
library('dplyr')
library('tibble')
library('readr')

registerDoMC(8)

load('saoGO.RData')
load('calKEGG.RData')
res <- read_csv('CAL_DEG_whole_k.csv')

## remove 0 terms
saoGO %<>% `[`(sapply(saoGO, length) > 0)
saoKEGG %<>% `[`(sapply(saoKEGG, length) > 0)

## padj < 0.05 & |log2FC| > 1
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
}
GOMat <- as.data.frame(GOMat)
GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE)
GOTestWithCat <- GOTestWithCat[!is.na(GOTestWithCat$ontology), ]

## add ablog2FC
goSub <- saoGO[match(GOTestWithCat[, 1], names(saoGO))]
abLogFC <- sapply(goSub, function(x) {
  eachFC <- res[match(x, res$GeneID), 'log2FoldChange']
  return(mean(abs(eachFC), na.rm = TRUE))
})
GOTestWithCat$abLogFC <- abLogFC

## deal with NA and select BP MF and CC
termCat <- c('BP', 'MF', 'CC')
for (i in termCat) {
  write.csv(GOTestWithCat[GOTestWithCat$ontology == i, ],
            paste0('degseq4h_FC2_', i, '_withcat.csv'))
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
################################################################

getwd()
setwd("C:/Users/bretl/Desktop/NDICG_practical_8/")

#BiocManager::install("airway")
browseVignettes("airway")
library(airway)
data("airway")
DS.airway <- airway

library("tidyverse")
library("tximport")
library("DESeq2")
dir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(dir)
list.files(file.path(dir, "quants"))

csvfile <- file.path(dir, "sample_table.csv")
coldata <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
coldata

coldata <- coldata[1:2,]
coldata$names <- coldata$Run
coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf.gz")
file.exists(coldata$files)

BiocManager::install("tximeta")
library("tximeta")
se <- tximeta(coldata)
dim(se)
head(rownames(se))
gse <- summarizeToGene(se)
data(gse)
gse
assayNames(gse)
head(assay(gse), 3)
colSums(assay(gse))
rowRanges(gse)
seqinfo(rowRanges(gse))
colData(gse)
gse$donor
gse$condition
gse$cell <- gse$donor
gse$dex <- gse$condition
levels(gse$dex)
levels(gse$dex) <- c("untrt", "trt")
library("magrittr")
gse$dex %<>% relevel("untrt")
gse$dex
gse$dex <- relevel(gse$dex, "untrt")
round( colSums(assay(gse)) / 1e6, 1 )
dds <- DESeqDataSet(gse, design = ~ cell + dex)
countdata <- round(assays(gse)[["counts"]])
head(countdata, 3)
coldata <- colData(gse)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ cell + dex)
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)
keep <- rowSums(counts(dds) >= 10) >= 3
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
BiocManager::install("vsn")
BiocManager::install("hexbin")
library("vsn")
library("hexbin")
meanSdPlot(cts, ranks = FALSE)
log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
library("dplyr")
library("ggplot2")
dds <- estimateSizeFactors(dds)
df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


sampleDists <- dist(t(assay(vsd)))
sampleDists
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

BiocManager::install("PoiClaClu")
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$dex, dds$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)
plotPCA(vsd, intgroup = c("dex", "cell"))
pcaData <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
BiocManager::install("glmpca")
library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$dex <- dds$dex
gpca.dat$cell <- dds$cell
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = dex, shape = cell)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")
dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds, contrast=c("dex","trt","untrt"))
mcols(res, use.names = TRUE)
summary(res)
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)
results(dds, contrast = c("cell", "N061011", "N61311"))
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(res$padj < 0.1, na.rm=TRUE)
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("dex"))
BiocManager::install("ggbeeswarm")
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex","cell"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
BiocManager::install("apeglm")
library("apeglm")
resultsNames(dds)
res <- lfcShrink(dds, coef="dex_trt_vs_untrt", type="apeglm")
plotMA(res, ylim = c(-5, 5))
res.noshr <- results(dds, name="dex_trt_vs_untrt")
plotMA(res.noshr, ylim = c(-5, 5))
plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
BiocManager::install("genefilter")
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
pheatmap(mat, annotation_col = anno)
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
  mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
        ylab = "fraction of small p values")

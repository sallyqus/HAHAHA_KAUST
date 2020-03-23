## EdgeR for DEG analysis
##Load the following R libraries in your script
library(limma)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(reshape2)
library(GGally)
library(EnhancedVolcano)
library(GO.db)
library(gprofiler2)
library(clusterProfiler)

## set your own directory as working directory
setwd("/Volumes/GoogleDrive/My Drive/Course_KAUST/B 390M Computational BioSci/2020Lec/assign2/B322/BulkRNASeqAnalysis/")

mycounts<-read.csv("GSE50760_counts_gene.csv",header = TRUE)
metadata <- read.table(file = 'GSE50760_meta.tsv', sep = '\t', header = TRUE)

rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
head(mycounts)

## Map from Ensembl gene id to gene symbol
ensg <- sub("\\..*", "", rownames(mycounts)) # remove version number
sym <- mapIds(org.Hs.eg.db, keys=ensg, column="SYMBOL", keytype="ENSEMBL")
gene <- data.frame(ENSGID=ensg, SYMBOL=sym, stringsAsFactors=F)
dim(gene)

t_ <- factor(metadata$Tissue, levels = c("primary_colorectal_cancer","normal_colon", "metastasized_cancer_to_liver"))
design <- model.matrix(~t_)
           
# Question 1 :
#   Explore the data
# 1a. remove genes with no expression for all samples
# 1b. boxplot of library size for Tissue group and Subjects
# 1c. filter genes by expression before normalization
# 1d. Look at the difference in cpm of raw expression and filtered expression values by density plots

# 1a. remove genes with no expression for all samples
y_raw <- DGEList(counts = mycounts, genes = gene) 
rownames(design) <- colnames(y_raw)
dim(y_raw)

na <- is.na(y_raw$genes$SYMBOL)
y <- y_raw[!na,]
dim(y)

d <- duplicated(y$genes$SYMBOL)
y <- y[!d,]
dim(y)

keep <- rowSums(y$counts) >= 1
y <- y[keep,]  
dim(y)

keep <- rowSums(cpm(y)>1) >= 1;
y <- y[keep,]  
dim(y)

# 1b. boxplot of library size for Tissue group and Subjects
y$samples$lib.size <- colSums(y$counts);
barplot(y$samples$lib.size*1e-6, names=1:54, ylab="Library size (millions)")


# 1c. Here, a gene is only retained if it is expressed at a minimum level:
# The filterByExpr function in the edgeR package provides an automatic way to filter genes, 
# while keeping as many genes as possible with worthwhile counts.
keep <- filterByExpr(y, group = t_)
summary(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)

# 1d. Look at the difference in cpm of raw expression and filtered expression values by density plots
L <- mean(y$samples$lib.size) * 1e-6
M <- median(y$samples$lib.size) * 1e-6
c(L, M)
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(y_raw)
lcpm_raw <- cpm(y_raw, log=TRUE)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm_raw[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm_raw[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
# legend("topright", colnames(cpm(y)), text.col=col, bty="n")
lcpm <- cpm(y, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
# legend("topright", colnames(cpm(y)), text.col=col, bty="n")

# Question 2 :
#   Data normalisation by TMM method
y <- calcNormFactors(y, method = "TMM");
y$samples
y$samples$norm.factors

par(mfrow=c(1,2))
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")

lcpm <- cpm(y, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

# Question 3 :
#   Run multidimensionality reduction like PCA or MDS
# Find out the outlier or not properly grouped samples and remove them for downstream analysis
par(mfrow=c(1,2))
plotMDS(y, col=rep(1:3, each=18));
title(main="Sample Name")

col.group <- t_
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=t_, col=col.group)
title(main="Tissue Group")


#remove SRR975588
y$counts <- y$counts[ ,-c(38)]
y$samples <- y$samples[-c(38), ]
design <- design[-38, ]

# Question 4 :
#   Differential expression for two contrasts a. PrimaryColon Vs Normal, b. MetastasisColon Vs PrimaryColon
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y);

fit <- glmFit(y, design, robust=TRUE)
## a. PrimaryColon Vs Normal
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt);
top

## b. MetastasisColon Vs PrimaryColon
lrt_2 <- glmLRT(fit, coef=3)
top_2 <- topTags(lrt_2);
top_2

# Question 5 :
#   Find differentially expressed genes at p value < 0.01 NS lfc=log2(4)
summary(de <- decideTestsDGE(lrt,p.value=0.01, lfc = log2(4) )); # default FDR = 0.05
detags <- rownames(y)[as.logical(de)];
plotSmear(lrt, de.tags=detags);
abline(h=c(-1, 1), col="blue");

summary(de_2 <- decideTestsDGE(lrt_2,p.value=0.01, lfc = log2(4) )); # default FDR = 0.05
detags_2 <- rownames(y)[as.logical(de_2)];
plotSmear(lrt_2, de.tags=detags_2);
abline(h=c(-1, 1), col="blue");


library(gplots)
i <- which(lrt$genes$ENSGID %in% top$table$ENSGID)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=lrt$genes$SYMBOL[i], labCol=t_, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

ma <- lcpm[which(lrt$genes$ENSGID %in% top$table$ENSGID),]

i_2 <- which(lrt_2$genes$ENSGID %in% top_2$table$ENSGID)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i_2,], scale="row",
          labRow=lrt_2$genes$SYMBOL[i_2], labCol=t_, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")


ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


# Question 6 :
#   Find the functional analysis of top 10 differentially expressed genes in case a and case b 

pdf(file=paste("D3_up", "Genecount_plot.pdf", sep="_"),width = 20, height =20)
#layout <- matrix(1:15, nrow = 3, ncol=5, byrow = TRUE)
for (i in 1:31){
  d <- plotCounts(ddss, gene=l1[i], intgroup="condition", 
                  returnData=TRUE)
  assign(paste("p",i,sep=""),ggplot(d, aes(x=condition, y=count, color = condition)) + 
           geom_point(position=position_jitter(w=0.1,h=0)) + 
           scale_y_log10(breaks=c(25,100,400)) +
           ylab("counts") +
           ggtitle(l[i]))
}
layout <- matrix(1:15, nrow = 5, ncol=3, byrow = TRUE)
print(multiplot(plotlist = list(p1,p2,p3,p4,p5,p6,p7,p8,p8,p9,p10,p11,p12,p13,p14,p15), layout = layout))
dev.off()

o <- order(lrt$table$PValue)
#na <- is.na(lrt[o,]$genes$SYMBOL)
#lrt_tmp_1 <- lrt[!na,]
de_10_genes <- lrt[o[1:10],]$genes$SYMBOL
egSYMBOL <- toTable(org.Hs.egSYMBOL)
m <- match(de_10_genes, egSYMBOL$symbol)
de_10_genes <- egSYMBOL$gene_id[m]
ego <- enrichGO(
  gene = de_10_genes,
  OrgDb = 'org.Hs.eg.db',
  ont = "all",
  readable = TRUE
)
head(ego, 10)
barplot(ego, showCategory = 5, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")


o_2 <- order(lrt_2$table$PValue)
de_10_genes_2 <- lrt_2[o_2[1:10],]$genes$SYMBOL
m <- match(de_10_genes_2, egSYMBOL$symbol)
de_10_genes_2 <- egSYMBOL$gene_id[m]
ego_2 <- enrichGO(
  gene = de_10_genes_2,
  OrgDb = 'org.Hs.eg.db',
  ont = "all",
  readable = TRUE
)
head(ego_2)
barplot(ego_2, showCategory = 5, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

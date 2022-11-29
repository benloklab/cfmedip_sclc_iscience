#################################################################
# This script closely examines the PRIME-filtered SCLC cfDNA data
# It then performs consensus clustering on the PBL-filtered cfDNA SCLC data to identify clusters of interest
# Written by Sami Ul haq
# Nov 8, 2021

library(ggplot2)
library(DESeq2)
library(ConsensusClusterPlus)
library(reshape2)
library(AnnotationHub)
source("../data/GSEA_KEGG_Analysis_Functions.R")

load("../data/PCA.cfDNA.PRIME.filtered.RData")
pc1.variance <- "PC1: 57% variance"
pc2.variance <- "PC2: 4% variance"

#################################################################
# Figure 3A - Consensus clustering on PRIME filtered SCLC cfDNA data

# Consensus Clustering done on RPKM-transformed PRIME-filtered cfDNA data
load("../data/PRIME.filtered.rpkm.matrix.RData")

# log transformed to normalize data
PRIME.filtered.rpkm.matrix <- PRIME.filtered.rpkm.matrix + 0.01
PRIME.filtered.rpkm.matrix <- log(PRIME.filtered.rpkm.matrix)
# hist(PRIME.filtered.rpkm.matrix)

# Consensus clustering is done
mads <- apply(PRIME.filtered.rpkm.matrix, 1, mad)
# Top 5000 most variant regions choosen
PRIME.filtered.rpkm.matrix <- PRIME.filtered.rpkm.matrix[rev(order(mads))[1:5000],]
PRIME.filtered.rpkm.matrix = sweep(PRIME.filtered.rpkm.matrix,1, apply(PRIME.filtered.rpkm.matrix,1,median,na.rm=T))

title <- "PRIME filter"
results <- ConsensusClusterPlus(PRIME.filtered.rpkm.matrix, maxK=8, reps=1000, pItem=0.8, pFeature=1, title=title, clusterAlg="hc", distance="pearson", seed=1262118388.71279, plot="pdf")
# The cluster info corresponding to patient is saved
consensus.clusters.df <- data.frame(k2=factor(results[[2]]$consensusClass),
                                    k3=factor(results[[3]]$consensusClass),
                                    k4=factor(results[[4]]$consensusClass))

# the clustering results are saved
#save(consensus.clusters.df, file="../results/Consensus.Clusters.PRIME_filtered.rpkm.RData")

# The clustering info is added to the pca data
pca.data$clusters.2 <- consensus.clusters.df$k2
pca.data$clusters.3 <- consensus.clusters.df$k3
pca.data$clusters.4 <- consensus.clusters.df$k4

levels(pca.data$clusters.2) <- c("A", "B")
named.color <- c("A"="#274FD9", "B"="#E6BA40")


ggplot(data = pca.data, aes(PC1, PC2, color=clusters.2, shape=clusters.2) ) + 
  geom_point(size=6) +
  theme_bw() +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        title = element_text(size = 24),
        legend.text = element_text(size=18)) +
  xlab(pc1.variance) + 
  ylab(pc2.variance) +
  scale_color_manual(values=named.color)

ggsave(filename = "../results/Figure_3A.pdf", units="in", width=12, height=6)


#################################################################
# Figure 3B - DMR analysis of the consensus clusters

# Examining Consensus Clusters in SCLC cfDNA
# Consensus clustering was done in the PRIME-filtered SCLC cfDNA data. This identified 2 clusters in the data.

load("../data/PRIME.filtered.counts.matrix.RData")

clusters.2 <- consensus.clusters.df$k2
clusters.2 <- as.character(clusters.2)
clusters.2[clusters.2 %in% "1"] <- "A"
clusters.2[clusters.2 %in% "2"] <- "B"
annotation.df <- data.frame(condition=c( clusters.2 ))

# DMR analysis between Cluster 1 & Cluster 2 is done
dds <- DESeqDataSetFromMatrix(countData = PRIME.filtered.counts.matrix, annotation.df, ~condition)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
# plotPCA(vsd, intgroup=c("condition"), n=500)
PCA.consensus.cluster.DMR.analysis <- plotPCA(vsd, intgroup=c("condition"), returnData = TRUE )
resultsNames(dds)

cluster.B.vs.cluster.A <- results(dds)
# plotMA(cluster.B.vs.cluster.A, ylim=c(-4,4))

# Shrinkage of effect size (better for minimizing variability)
cluster.B.vs.cluster.A <- lfcShrink(dds, coef="condition_B_vs_A")
# plotMA(cluster.B.vs.cluster.A, ylim=c(-4,4))
# save(cluster.B.vs.cluster.A, file="DMR.consensus.cluster.B.vs.cluster.A.RData")

# The DESEQ2 analysis examines Consensus Cluster B vs A. This means log2(Cluster B / Cluster A)
# So, (+)Log2FC (and p-val < 0.05) = windows enriched in Cluster B
# and (-)Log2FC (and p-val < 0.05) = windows enriched in Cluster A

p.val.cutoff <- 0.05
cluster.a.fc.cutoff <- 2
cluster.b.fc.cutoff <- 1

# volcano plot is made for the DMR analysis
dsq2 <- data.frame(log2FC=cluster.B.vs.cluster.A$log2FoldChange, adjust.pval=cluster.B.vs.cluster.A$padj)
dsq2 <- na.omit(dsq2)

dsq2$col <- rep("NS", length(dsq2$log2FC))
dsq2$col[which(dsq2$adjust.pval < 0.05 & dsq2$log2FC > 2)] <- "Cluster B"
dsq2$col[which(dsq2$adjust.pval < 0.05 & dsq2$log2FC < -1)] <- "Cluster A"

named.color.vector <- c("NS" = "#BAC4C8", "Cluster A" = "#435AA9", "Cluster B" = "#E6BA40")

vulcan <- ggplot(dsq2, aes(x=log2FC, y=-log10(adjust.pval), color=col)) +
  geom_point(alpha=0.5, size=5) +
  scale_color_manual(values = named.color.vector) +
  theme_bw() +
  geom_vline(xintercept= cluster.a.fc.cutoff, linetype=1, size=2, color="#B79896") +
  geom_vline(xintercept= -cluster.b.fc.cutoff, linetype=1, size=2, color="#B79896") +
  geom_hline(yintercept = -log10(0.05), linetype=1, size=2, color="#B79896") +
  theme(axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        legend.position = "bottom")

ggsave(plot=vulcan, filename = "../results/Figure_3B.tiff", units = "in", width = 8, height = 7, dpi = 300)
# the plot is deleted after saving because it is memory intensive
rm(vulcan, dsq2)

#######################################################################
# Figure 3C - KEGG Pathway analysis

# The following examines the characteristics of the windows that are hypermethylated in Cluster B vs A. Specifically, to get an idea of genomic locations commonly methylated (i.e. promoter regions, gene bodies, etc.)

regions.in.clusterA <- rownames(cluster.B.vs.cluster.A)[which(cluster.B.vs.cluster.A$padj < p.val.cutoff & cluster.B.vs.cluster.A$log2FoldChange < -1)]
regions.in.clusterA.grange <- medip.regions.to.GRange.maker(regions.in.clusterA)

regions.in.clusterB <- rownames(cluster.B.vs.cluster.A)[which(cluster.B.vs.cluster.A$padj < p.val.cutoff & cluster.B.vs.cluster.A$log2FoldChange > 2)]
regions.in.clusterB.grange <- medip.regions.to.GRange.maker(regions.in.clusterB)

# Pathway/KEGG analysis of the methylated windows are examined in cluster A vs B
clusterA.genes <- medip.Grange.to.gene.symbols(regions.in.clusterA.grange)
# cat("\nThere are ", length(clusterA.genes), "in cluster A\n")
clusterB.genes <- medip.Grange.to.gene.symbols(regions.in.clusterB.grange)
# cat("\nThere are ", length(clusterB.genes), "in cluster B\n")

clusterA.genes.kegg <- kegg.analysis(clusterA.genes)
barplot(clusterA.genes.kegg, showCategory=20)
ggsave(filename="../results/Figure_3C_cluster.A.pdf")
# blah <- setReadable(clusterA.genes.kegg, OrgDb = org.Hs.eg.db, keyType="UNIPROT")
# write.csv(blah, file="cluster.A.KEGG_genes.csv")


clusterB.genes.kegg <- kegg.analysis(clusterB.genes)
barplot(clusterB.genes.kegg, showCategory=20)
ggsave(filename="../results/Figure_3C_cluster.B.pdf")
# blah <- setReadable(clusterB.genes.kegg, OrgDb = org.Hs.eg.db, keyType="UNIPROT")
# write.csv(blah, file="cluster.B.KEGG_genes.csv")

################################################################################
# Figure 4A - Differences in CpG features methylation is examined 

annots = c("hg19_cpgs")
annotations = build_annotations(genome = 'hg19', annotations = annots)

regions.in.clusterA.annotated <- genomic.window.annotator(queryGrange = regions.in.clusterA.grange, subjectGrange = annotations, removeDuplicateRegions = FALSE)
# cpg features in cluster A are saved
table.of.clust.a.cpg.feat <- c(table(regions.in.clusterA.annotated$type))

regions.in.clusterB.annotated <- genomic.window.annotator(queryGrange = regions.in.clusterB.grange, subjectGrange = annotations, removeDuplicateRegions = FALSE)
# cpg features in cluster B are saved
table.of.clust.b.cpg.feat <- c(table(regions.in.clusterB.annotated$type))

features.clusterA <- table.of.clust.a.cpg.feat / sum(table.of.clust.a.cpg.feat)
features.clusterB <- table.of.clust.b.cpg.feat / sum(table.of.clust.b.cpg.feat)

# empty dataframe is made
temp.annots.holder <- data.frame(row.names = 1:length(features.clusterB), stringsAsFactors = FALSE)
temp.annots.holder$Features <- as.vector( names(features.clusterB))
temp.annots.holder$Cluster.A <- rep(0, length(features.clusterB) )
temp.annots.holder$Cluster.B <- rep(0, length(features.clusterB) )

#dataframe is populated with the appropriate frequencies of each annotation 
temp.annots.holder$Cluster.A[match( names(features.clusterA),  temp.annots.holder$Features)] <- features.clusterA
temp.annots.holder$Cluster.B[match( names(features.clusterB),  temp.annots.holder$Features)] <- features.clusterB

reform.annots.df <- melt(temp.annots.holder,
                         measure.vars = c("Cluster.A", "Cluster.B"), 
                         variable.name = "Clusters")

cpg.graph <- ggplot(reform.annots.df, aes(Clusters, value, fill=as.factor(Features))) + 
  geom_col(width = 0.8) +
  theme_classic() +
  labs(fill="CpG Feature") +
  theme(axis.text.x = element_text()) +
  scale_fill_manual(values = c("hg19_cpg_inter"="#66c2a5", 
                               "hg19_cpg_islands"="#fc8d62", 
                               "hg19_cpg_shelves"="#8da0cb", 
                               "hg19_cpg_shores"="#e78ac3")) +
  scale_x_discrete(labels=c("Cluster A", "Cluster B")) +
  ylab("Proportion of Significant DMRs") +
  theme(legend.position = "right",
        axis.text.x = element_text(size=18, angle = 0),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        title = element_text(size = 24),
        legend.text = element_text(size=18)) +
  scale_y_continuous(labels = scales::percent)


ggsave(plot=cpg.graph, filename="../results/Figure_4A.pdf", width=10, height=8, units = "in")


##############################################################################
# Figure 4B - Differences in methylation of repeats regions is examined

load("../data/ucsc.repeatmasker.annotations.RData")
features.to.consider <- c("LINE", "SINE", "LTR", "Satellite", "Retroposon")
ucsc.repeatmasker.annotations <- subset(ucsc.repeatmasker.annotations, ucsc.repeatmasker.annotations$repeat.class %in% features.to.consider)

regions.in.clusterA.annotated <- genomic.window.annotator(regions.in.clusterA.grange, ucsc.repeatmasker.annotations)
table.of.clust.a.repeats <- c(table(regions.in.clusterA.annotated$repeat.class))
len.of.table <- length(table.of.clust.a.repeats) + 1
table.of.clust.a.repeats[len.of.table] <- length(regions.in.clusterA.grange) -
  sum(table(regions.in.clusterA.annotated$repeat.class))
names(table.of.clust.a.repeats)[ len.of.table ] <- "A Non-repeats feature"

regions.in.clusterB.annotated <- genomic.window.annotator(regions.in.clusterB.grange, ucsc.repeatmasker.annotations)
table.of.clust.b.repeats <- c(table(regions.in.clusterB.annotated$repeat.class))
len.of.table <- length(table.of.clust.b.repeats) + 1
table.of.clust.b.repeats[len.of.table] <- length(regions.in.clusterB.grange) -
  sum(table(regions.in.clusterB.annotated$repeat.class))
names(table.of.clust.b.repeats)[ len.of.table ] <- "A Non-repeats feature"


features.clusterA <- table.of.clust.a.repeats / length(regions.in.clusterA.grange)
features.clusterB <- table.of.clust.b.repeats / length(regions.in.clusterB.grange)


# empty dataframe is made
temp.annots.holder <- data.frame(row.names = 1:length(features.clusterB), stringsAsFactors = FALSE)
temp.annots.holder$Features <- as.vector( names(features.clusterB))
temp.annots.holder$Cluster.A <- rep(0, length(features.clusterB) )
temp.annots.holder$Cluster.B <- rep(0, length(features.clusterB) )

#dataframe is populated with the appropriate frequencies of each annotation 
temp.annots.holder$Cluster.A[match( names(features.clusterA),  temp.annots.holder$Features)] <- features.clusterA
temp.annots.holder$Cluster.B[match( names(features.clusterB),  temp.annots.holder$Features)] <- features.clusterB

# change order of data rows
temp.annots.holder <- temp.annots.holder[c(6,1,2,3,4,5) ,]

reform.annots.df <- melt(temp.annots.holder,
                         measure.vars = c("Cluster.A", "Cluster.B"), 
                         variable.name = "Clusters")


named.color.vector <- c("LINE"="#D460DB",
                        "LTR"="#85C1DE", 
                        "Retroposon"="#DEC6A8", 
                        "Satellite"="#DB9239",
                        "SINE"="#A9E080",
                        "A Non-repeats feature"="#e6e3e3")

repeats.feat <- ggplot(reform.annots.df, aes(Clusters, value, fill=as.factor(Features))) + 
  geom_col(width = 0.8) +
  theme_classic() +
  labs(fill="Repeats Feature") +
  theme() +
  ylab("Proportion of Significant DMRs") +
  theme(legend.position = "right",
        axis.text.x = element_text(size=18, angle = 0),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        title = element_text(size = 24),
        legend.text = element_text(size=18),
        legend.title = element_blank()) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = named.color.vector) +
  scale_x_discrete(labels=c("Cluster A", "Cluster B"))
  
  
ggsave(plot=repeats.feat, filename="../results/Figure_4B.pdf", width=10, height=8, units = "in")


###########################################################################
# Figure 4C - Differences in methylation of lncRNA regions is examined

load("../data/non.coding.gene.features.RData")
features.to.consider <- c("lnc_RNA")

# only HAVANA classified lncRNA are considered 
non.coding.gene.features <- subset(non.coding.gene.features, non.coding.gene.features$source %in% "havana")
non.coding.gene.features <- subset(non.coding.gene.features, non.coding.gene.features$biotype %in% "lncRNA")
non.coding.gene.features <- subset(non.coding.gene.features, non.coding.gene.features$type %in% features.to.consider)


regions.in.clusterA.annotated <- genomic.window.annotator(regions.in.clusterA.grange, non.coding.gene.features, removeDuplicateRegions = FALSE)
regions.in.clusterB.annotated <- genomic.window.annotator(regions.in.clusterB.grange, non.coding.gene.features, removeDuplicateRegions = FALSE)


lncRNA.cluster.A <- table(regions.in.clusterA.annotated$type)[4] / length(regions.in.clusterA.grange)
non.lncRNA.windows.cluster.A <- 1 - lncRNA.cluster.A
lncRNA.cluster.B <- table(regions.in.clusterB.annotated$type)[4] / length(regions.in.clusterB.grange)
non.lncRNA.windows.cluster.B <- 1 - lncRNA.cluster.B


lncrna.by.cluster <- data.frame(Features = c("lncRNA", "lncRNA", "A Non-lncRNA", "A Non-lncRNA"),
                                Clusters = c("Cluster A", "Cluster B", "Cluster A", "Cluster B"),
                                value = c(lncRNA.cluster.A, lncRNA.cluster.B, non.lncRNA.windows.cluster.A, non.lncRNA.windows.cluster.B))


named.color.vector <- c("lncRNA"="#fc0388", "A Non-lncRNA"="#e6e3e3")

# lncRNA features
lncRNA.feat <- ggplot(lncrna.by.cluster, aes(Clusters, value, fill=as.factor(Features))) + 
  geom_col(width = 0.8) +
  theme_classic() +
  labs(fill="cluster") +
  theme(axis.text.x = element_text()) +
  scale_fill_manual(values = named.color.vector) +
  ylab("Proportion of DMRs") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size=18, angle = 0),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        title = element_text(size = 24),
        legend.text = element_text(size=18),
        legend.title = element_blank()) +
  scale_y_continuous(labels = scales::percent)

ggsave(plot=lncRNA.feat, filename="../results/Figure_4C.pdf", width=7, height=8, units = "in")

#################################################################
# This script closely correlates the PRIME-filtered cfDNA data to
# various features of interest in the genome
# Written by Sami Ul haq
# April  2022

library(ggplot2)
library(DESeq2)
library(annotatr)
source("../data/GSEA_KEGG_Analysis_Functions.R")


# The principal component anlayses in the later section of the code were used to generate the PC1 variance values
# correlating PC1 variance vs number of windows used
corr.pc.vs.windows <- data.frame(features = c("Promoters", "TFs", "miRNA", "lncRNA", "SINEs", "Satellites", "SNVs", "Exons", "CTCF-site", "LTRs", "LINEs"),
                                 num.windows = c(11702, 1626, 65, 80236, 64027, 1000, 919, 34042, 10715, 25917, 59901),
                                 PC1.var = c(0.45, 0.35, 0.19, 0.58, 0.58, 0.16, 0.26, 0.57, 0.49, 0.54, 0.58))

named.color <- c("Promoters" = "#a6cee3", "TFs" = "#b2df8a", 
                 "miRNA" = "#33a02c", "lncRNA" = "#fb9f9d", 
                 "SINEs" = "#e31a1c", "Satellites" = "#fdbc68", 
                 "SNVs" = "#ff7f00", "Exons" = "#cab2d6", 
                 "CTCF-site" = "#6a3d9a", "LTRs" = "#e6e600", "LINEs" = "#d0672f")

ggplot(corr.pc.vs.windows, aes(num.windows, PC1.var, color=features)) +
  geom_point(size=6) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        title = element_text(size = 24),
        legend.text = element_text(size=18)) +
  xlab("Number of PRIME-filtered windows") +
  ylab("PC1 Variance (%)") +
  scale_y_continuous(label=scales::percent) +
  scale_color_manual(values = named.color)
  
ggsave(filename="../results/Figure_4E.pdf", units="in",width = 10, height=6)



########################################################################
# Figure 4F & G - PCA of lncRNA, LINEs and SINEs
# Annotation data is read in

# Regulatory features
load("../data/regulatory.features.RData")

regulatory.features.ensembl.promoter <- subset(regulatory.features.ensembl, regulatory.features.ensembl$feature_type %in% "Promoter")
regulatory.features.ensembl.CTCF <- subset(regulatory.features.ensembl, regulatory.features.ensembl$feature_type %in% "CTCF Binding Site")
regulatory.features.ensembl.TF <- subset(regulatory.features.ensembl, regulatory.features.ensembl$feature_type %in% "TF binding site")

# Clinically relevant variants
load("../data/clinically.relevant.variants.RData")
clin.relevant.variants.SNVs <- subset(clin.relevant.variants, clin.relevant.variants$type %in% "SNV")

# Non coding regions
load("../data/non.coding.gene.features.RData")
non.coding.lnc.RNA <- subset(non.coding.gene.features, non.coding.gene.features$type %in% "lnc_RNA")
non.coding.miRNA <- subset(non.coding.gene.features, non.coding.gene.features$type %in% "miRNA")

# Repeats regions
load("../data/ucsc.repeatmasker.annotations.RData")
LINEs.only <- subset(ucsc.repeatmasker.annotations, ucsc.repeatmasker.annotations$repeat.class %in% "LINE")
SINEs.only <- subset(ucsc.repeatmasker.annotations, ucsc.repeatmasker.annotations$repeat.class %in% "SINE")
LTRs.only <- subset(ucsc.repeatmasker.annotations, ucsc.repeatmasker.annotations$repeat.class %in% "LTR")
Satellites.only <- subset(ucsc.repeatmasker.annotations, ucsc.repeatmasker.annotations$repeat.class %in% "Satellite")

# exons for protein coding genes
annots = c("hg19_genes_exons")
exons.only = build_annotations(genome = 'hg19', annotations = annots)


###################################################################################
features.pca.function <- function(sample.conditions, counts.matrix, file.to.name, special.case=FALSE) {
  #' 
  #' @description
  #' This function takes as input a counts matrix subsetted to specific features
  #' and performs PCA using DESEQ2's native pca function. The function
  #' returns a PCA plot (PDF & RData)
  #' 
  #' few rows with mean normalized count > 5 then use special.case = TRUE
  #' 
  #' 
  
  annotation.df <- data.frame(condition=c( sample.conditions ))
  
  dds <- DESeqDataSetFromMatrix(countData = counts.matrix, annotation.df, ~condition)
  dds <- DESeq(dds)
  
  # few rows with mean normalized count > 5
  if(special.case) {
    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  } else {
    vsd <- vst(dds, blind=FALSE)
  }
  
  plot.to.return <- plotPCA(vsd, intgroup=c("condition"))
  
  return(plot.to.return)
  
}

pca.plotter.and.saver <- function(pca.plot.object, file.to.name) {
  #' 
  #' @description 
  #' This function takes as input a ggplot PCA plot object and filename, plots it
  #' and saves the pca data using the filename provided
  #' 
  
  pca.data <- pca.plot.object$data
  save(pca.data, file=paste0(file.to.name, ".RData"))
  ggsave(plot = pca.plot.object, filename=paste0(file.to.name, ".pdf"))
  
}

# SCLC cfDNA data using PRIME filter
load("../data/PRIME.filtered.counts.matrix.RData")
#dim(PRIME.filtered.counts.matrix)

# the windows in the matrix are converted to a GRange
matrix.windows.grange <- medip.regions.to.GRange.maker(rownames(PRIME.filtered.counts.matrix))

load("../data/consensus.clusters.prime.filtered.RData")
sample.conditions <- consensus.clusters.df$k2


##############################################################
# Examining the overall cfDNA data after PRIME filter 

#file.to.name <- "PRIME.filtered"
#PRIME.filtered.counts <- features.pca.function(sample.conditions, PRIME.filtered.#counts.matrix, file.to.name)
#pca.plotter.and.saver(PRIME.filtered.counts , file.to.name)


# Examining the data after PRIME filter for lnc RNA
# 80236 windows
# finds the windows that overlap with lnc RNA only
lnc.RNA.only.windows <- windows.overlapping.with.features(matrix.windows.grange, non.coding.lnc.RNA)
matri.lnc.RNA <- PRIME.filtered.counts.matrix[which(rownames(PRIME.filtered.counts.matrix) %in% lnc.RNA.only.windows), ]
dim(matri.lnc.RNA)

file.to.name <- "../results/Figure_4F"
PRIME.filtered.counts.lncRNA <- features.pca.function(sample.conditions, matri.lnc.RNA, file.to.name)
pca.plotter.and.saver(PRIME.filtered.counts.lncRNA , file.to.name)



# Examining the data after PRIME filter for SINEs
# 64027 windows
# finds the windows that overlap with SINEs only
SINEs.only.windows <- windows.overlapping.with.features(matrix.windows.grange, SINEs.only)
matri.SINEs.only <- PRIME.filtered.counts.matrix[which(rownames(PRIME.filtered.counts.matrix) %in% SINEs.only.windows), ]
dim(matri.SINEs.only)

file.to.name <- "../results/Figure_4G"
PRIME.filtered.counts.SINEs <- features.pca.function(sample.conditions, matri.SINEs.only, file.to.name)
pca.plotter.and.saver(PRIME.filtered.counts.SINEs , file.to.name)


# Examining the data after PRIME filter for LINEs
# 59901 windows
# finds the windows that overlap with LINES only
LINEs.only.windows <- windows.overlapping.with.features(matrix.windows.grange, LINEs.only)
matri.LINEs.only <- PRIME.filtered.counts.matrix[which(rownames(PRIME.filtered.counts.matrix) %in% LINEs.only.windows), ]
dim(matri.LINEs.only)

file.to.name <- "../results/Figure_4H"
PRIME.filtered.counts.LINEs <- features.pca.function(sample.conditions, matri.LINEs.only, file.to.name)
pca.plotter.and.saver(PRIME.filtered.counts.LINEs , file.to.name)



######################################
# The following code examines other PRIME-filtered features of interest (promoters, exons, CTCFs, LTRs)

# Examining the data after PRIME filter for promoters
# 11702 windows
# finds the windows that are Promoters only & subsets the matrix
#promoters.only.windows <- windows.overlapping.with.features(matrix.windows.grange, regulatory.features.ensembl.promoter)
#matri.promoters <- PRIME.filtered.counts.matrix[which(rownames(PRIME.filtered.counts.matrix) %in% promoters.only.windows), ]
#dim(matri.promoters)

#file.to.name <- "PRIME.filtered.promoters"
#PRIME.filtered.counts.promoters <- features.pca.function(sample.conditions, matri.promoters, file.to.name)
#pca.plotter.and.saver(PRIME.filtered.counts.promoters , file.to.name)


# Examining the data after PRIME filter for exons
# 34042 windows
# finds the windows that are Promoters only & subsets the matrix
#exons.only.windows <- windows.overlapping.with.features(matrix.windows.grange, exons.only)
#matri.exons <- PRIME.filtered.counts.matrix[which(rownames(PRIME.filtered.counts.matrix) %in% exons.only.windows), ]
#dim(matri.exons)

#file.to.name <- "PRIME.filtered.exons"
#PRIME.filtered.counts.exons <- features.pca.function(sample.conditions, matri.exons, file.to.name)
#pca.plotter.and.saver(PRIME.filtered.counts.exons , file.to.name)


# Examining the data after PRIME filter for CTCF
# 10715 windows
# finds the windows that are CTCF (11-zinc finger protein or CCCTC-binding factor) only
#CTCF.only.windows <- windows.overlapping.with.features(matrix.windows.grange, regulatory.features.ensembl.CTCF)
#matri.CTCF <- PRIME.filtered.counts.matrix[which(rownames(PRIME.filtered.counts.matrix) %in% CTCF.only.windows), ]
#dim(matri.CTCF)

#file.to.name <- "PRIME.filtered.CTCF"
#PRIME.filtered.counts.CTCF <- features.pca.function(sample.conditions, matri.CTCF, file.to.name)
#pca.plotter.and.saver(PRIME.filtered.counts.CTCF , file.to.name)


# Examining the data after PRIME filter for TF
# 1626 windows
# finds the windows that are TF (transcription factors) only
#TF.only.windows <- windows.overlapping.with.features(matrix.windows.grange, regulatory.features.ensembl.TF)
#matri.TF <- PRIME.filtered.counts.matrix[which(rownames(PRIME.filtered.counts.matrix) %in% TF.only.windows), ]
#dim(matri.TF)

#file.to.name <- "PRIME.filtered.TF"
#PRIME.filtered.counts.TF <- features.pca.function(sample.conditions, matri.TF, file.to.name)
#save(PRIME.filtered.counts.TF, file=paste0(file.to.name, ".RData"))



# Examining the data after PRIME filter for SNVs
# 919 windows
# finds the windows that overlap with SNVs only
#SNVs.only.windows <- windows.overlapping.with.features(matrix.windows.grange, clin.relevant.variants.SNVs)
#matri.SNV <- PRIME.filtered.counts.matrix[which(rownames(PRIME.filtered.counts.matrix) %in% SNVs.only.windows), ]
#dim(matri.SNV)

#file.to.name <- "PRIME.filtered.SNVs"
#PRIME.filtered.counts.SNV <- features.pca.function(sample.conditions, matri.SNV, file.to.name, special.case = TRUE)
#pca.plotter.and.saver(PRIME.filtered.counts.SNV , file.to.name)


# Examining the data after PRIME filter for miRNA
# 65 windows
# finds the windows that overlap with miRNA only
#miRNA.RNA.only.windows <- windows.overlapping.with.features(matrix.windows.grange, non.coding.miRNA)
#matri.miRNA <- PRIME.filtered.counts.matrix[which(rownames(PRIME.filtered.counts.matrix) %in% miRNA.RNA.only.windows), ]
#dim(matri.miRNA)

#file.to.name <- "PRIME.filtered.miRNA"
#PRIME.filtered.counts.miRNA <- features.pca.function(sample.conditions, matri.miRNA, file.to.name, special.case = TRUE)
#pca.plotter.and.saver(PRIME.filtered.counts.miRNA , file.to.name)


# Examining the data after PRIME filter for Satellites
# 1000 windows
# finds the windows that overlap with Satellites only
#Satellites.only.windows <- windows.overlapping.with.features(matrix.windows.grange, Satellites.only)
#matri.Satellites.only <- PRIME.filtered.counts.matrix[which(rownames(PRIME.filtered.counts.matrix) %in% Satellites.only.windows), ]
#dim(matri.Satellites.only)

#file.to.name <- "PRIME.filtered.Satellites.only"
#PRIME.filtered.counts.Satellites <- features.pca.function(sample.conditions, matri.Satellites.only, file.to.name, special.case = TRUE)
#pca.plotter.and.saver(PRIME.filtered.counts.Satellites , file.to.name)


# Examining the data after PRIME filter for LTRs
# 25917 windows
# finds the windows that overlap with LTRs only
#LTRs.only.windows <- windows.overlapping.with.features(matrix.windows.grange, LTRs.only)
#matri.LTRs.only <- PRIME.filtered.counts.matrix[which(rownames(PRIME.filtered.counts.matrix) %in% LTRs.only.windows), ]
#dim(matri.LTRs.only)

#file.to.name <- "PRIME.filtered.LTRs.only"
#PRIME.filtered.counts.LTRs <- features.pca.function(sample.conditions, matri.LTRs.only, file.to.name)
#pca.plotter.and.saver(PRIME.filtered.counts.LTRs , file.to.name)

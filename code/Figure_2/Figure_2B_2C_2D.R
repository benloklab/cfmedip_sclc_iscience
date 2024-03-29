#########################################
# Figure 2B,C,D plots
# This script plots PCA for SCLC PBLs vs cfDNA,
# unfiltered SCLC cfDNA and 
# PRIME-filtered SCLC cfDNA PCA plots
#
# Written by Sami Ul Haq (Dec 2021)

library(ggplot2)

#######################################################
# Figure 2B - SCLC PBLs vs cfDNA

# The PCA was generated by examining 300bp ENCODE-blacklisted filtered, whole-genome methylation profiles between SCLC cfDNA and paired patient PBLs using limma on an HPC node with 180GB memory.
# See Zenodo link for downloading data:
# Ul Haq, Sami, Tsao, Ming S., Cabanero, Michael, de Carvalho, Daniel, Liu, Geoffrey, Bratman, Scott V., & Lok, Benjamin H. (2022). Processed counts data of cfMeDIP-seq profiles of small cell lung cancer patients [Data set]. In Cell Press iScience (Vol. 25, Number 12). Zenodo. https://doi.org/10.5281/zenodo.7235989

load("../data/PCA.plot.cfDNA.vs.pbl.RData")

# the labels are made more readable
sample.type <- substr(names(PCA.plot$x), 10, 14)

mds.df <- data.frame(PC.1 = PCA.plot$x,
                     PC.2 = PCA.plot$y,
                     sample = sample.type)

MDS.plot.cfDNA.vs.pbl <- ggplot(mds.df, aes(x=PC.1, y=PC.2, shape=sample, color=sample)) +
                                geom_point(size=4) +
                                theme_bw() +
                                theme(axis.text.x = element_text(size = 18),
                                axis.text.y = element_text(size = 18),
                                axis.title.y = element_text(size = 20),
                                axis.title.x = element_text(size = 20),
                                title = element_text(size = 24),
                                legend.text = element_text(size=18))


ggsave(MDS.plot.cfDNA.vs.pbl, file="../results/Figure_2B.pdf", width=10, height=8, units="in")


#######################################################
# Figure 2C - SCLC cfDNA (unfiltered)

# The PCA was generated by examining 300bp ENCODE-blacklisted filtered, whole-genome methylation profiles between SCLC cfDNA patients using DESeq2 on an HPC node with 180GB memory. The plotPCA() function in DESeq2 was used to generate the PCA and save the PCA data
# See Zenodo link for downloading data:
# Ul Haq, Sami, Tsao, Ming S., Cabanero, Michael, de Carvalho, Daniel, Liu, Geoffrey, Bratman, Scott V., & Lok, Benjamin H. (2022). Processed counts data of cfMeDIP-seq profiles of small cell lung cancer patients [Data set]. In Cell Press iScience (Vol. 25, Number 12). Zenodo. https://doi.org/10.5281/zenodo.7235989

load("../data/PCA.SCLC.cfDNA.unfiltered.RData")

pc1.variance <- "PC1: 40% variance"
pc2.variance <- "PC2: 7% variance"

ggplot(pca.data, aes(PC1, PC2)) +
  geom_point(size=3, color="red") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        title = element_text(size = 24),
        legend.text = element_text(size=18)) +
  xlab(pc1.variance) +
  ylab(pc2.variance)

ggsave(file="../results/Figure_2C.pdf", width=8, height=5, units="in")


#######################################################
# Figure 2D - PRIME-filtered SCLC cfDNA analysis

# The PCA was generated by examining 300bp ENCODE-blacklist filtered, PRIME-filtered cfDNA methylation profiles between SCLC cfDNA patients using DESeq2 on an HPC node with 60GB memory. The PRIME (Peripheral blood leukocyte methylation) windows were generated by generating MeDEStrand-infered beta values for whole-genome ENCODE-blacklist filtered methylation counts data for SCLC PBLs. Windows with a median beta value <= 0.2 across all PBL samples were kept. These windows were subsequently filtered for CG density >= 5 to account for 5meC antibody function. This was done on an HPC compute node with a memory of 180GB. These PRIME windows were examined in the SCLC cfDNA methylation profiles using DESeq2 on a compute node with 60GB memory and the plotPCA function was used.

# See Zenodo link for downloading data:
# Ul Haq, Sami, Tsao, Ming S., Cabanero, Michael, de Carvalho, Daniel, Liu, Geoffrey, Bratman, Scott V., & Lok, Benjamin H. (2022). Processed counts data of cfMeDIP-seq profiles of small cell lung cancer patients [Data set]. In Cell Press iScience (Vol. 25, Number 12). Zenodo. https://doi.org/10.5281/zenodo.7235989

load("../data/PCA.cfDNA.PRIME.filtered.RData")

pc1.variance <- "PC1: 57% variance"
pc2.variance <- "PC2: 4% variance"

ggplot(pca.data, aes(PC1, PC2)) +
  geom_point(size=3, color="red") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        title = element_text(size = 24),
        legend.text = element_text(size=18)) +
  xlab(pc1.variance) +
  ylab(pc2.variance)

ggsave(file="../results/Figure_2D.pdf", width=8, height=5, units="in")


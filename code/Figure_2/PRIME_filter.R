############################################################
# The following is the PRIME-filter cfDNA workflow
# Written by Sami Ul Haq
# Dec 2021
#
# This script is run on an HPC compute node with 180GB memory


library("MeDEStrand")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg19")
library("BiocGenerics")
library("S4Vectors")
library("IRanges")
library("GenomicRanges")
library("Biostrings")
library("XVector")
library("rtracklayer")
library("Rsamtools")
library("MEDIPS")
library("Repitools")

############################################################
# This part infers beta-values from PBL MeDIP counts data

BSgenome="BSgenome.Hsapiens.UCSC.hg19"
uniq = 1
extend = 200
shift = 0
ws = 300
chr.select = paste0('chr', c(1:22) )

# PBL BAM files
bam.files.names <- Sys.glob("path_to_bams")

# creates 300bp tiled genome for chromosomes 1 - 22
genome_300bp <- data.frame(genomeBlocks(BSgenome.Hsapiens.UCSC.hg19, chrs=seqnames(BSgenome.Hsapiens.UCSC.hg19)[1:22], width=300))
genome_300bp$name_format <- paste(genome_300bp$seqnames, genome_300bp$start, genome_300bp$end, sep=".")

# this matrix will hold all medestrand beta values
matrix.of.pbls <- matrix(nrow=length(genome_300bp$name_format))
# the rownames correspond to windows
rownames(matrix.of.pbls) <- genome_300bp$name_format

# each bam file is examined and the beta-value is infered using MeDEStrand and then added as a column to the holder matrix
for(each.bam.file in bam.files.names) {
  MeDIP_seq = MeDEStrand.createSet(file=each.bam.file, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws, chr.select=chr.select, paired = T)

  #  count CpG pattern in the bins
  CS = MeDEStrand.countCG(pattern="CG", refObj=MeDIP_seq)
  result.methylation = MeDEStrand.binMethyl(MSetInput = MeDIP_seq, CSet = CS, Granges = FALSE)

  # each sample is added
  matrix.of.pbls <- cbind(matrix.of.pbls, result.methylation)
}


############################################################
# this examines the median beta values per window
median.beta.vals <- apply(matrix.of.pbls, MARGIN=1, FUN = median)
names(median.beta.vals) <- rownames(matrix.of.pbls)

# Filter out blacklist windows
load("hg19.encode.blacklist.windows.RData")
# matrix is removed for blacklist windows
median.beta.vals <- median.beta.vals[ !(names(median.beta.vals) %in% hg19.encode.blacklist.windows) ]

# selects windows with beta values less than 0.2
PBL.depleted.windows <- median.beta.vals[which(median.beta.vals < 0.2)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.2.RData")

############################################################
# This script calculates the # of CGs per 300 bp window
# Written by Sami Ul Haq

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

running.total.cgs.per.window <- c()

# The Number of CGs per window is calculated for all windows
debug.counter <- 0
cat("\nCalculating CGs per window\n")

buncha.windows <- vector.containing.300bp.windows.in.genome

for (each.window in buncha.windows) {
  # the window name is stored
  window.name <- each.window
  # the chromosome, start, and end is stored
  window.coordinates <- unlist(strsplit(window.name, split="[.]"))
  # a dataframe followed by a GRange is made (convenient data structure for subsequent analysis)
  temp.df <- data.frame(seqnames=window.coordinates[1],
                        start=window.coordinates[2],
                        end=window.coordinates[3])
  temp.grange <- makeGRangesFromDataFrame(temp.df)
  rm(temp.df)

  # the fragment of DNA associated with the above window is pulled
  frag.of.DNA <- getSeq(BSgenome.Hsapiens.UCSC.hg19, temp.grange, as.character=TRUE)

  num.of.cgs <- as.numeric(vcountPattern("CG",frag.of.DNA))
  running.total.cgs.per.window <- c(running.total.cgs.per.window, num.of.cgs)

  debug.counter <- debug.counter + 1

  if(debug.counter %% 1000 == 0) {
    cat("\nCalculating window ", debug.counter, "\n")
  }
}
cat("\nFinished calculating CGs per window\n")

names(running.total.cgs.per.window) <- buncha.windows

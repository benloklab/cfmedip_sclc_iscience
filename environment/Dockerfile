# hash:sha256:2e51abcb043ac0725f88ea54ae23c4297fdf673ad0c396b55d027d3f6d1c433f
FROM registry.codeocean.com/codeocean/r-studio:2022.07.0-548-r4.2.1-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys \
        0xAD2323F17326AE31401037733E05EBFF05441C52 \
    && apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && add-apt-repository -y 'deb http://deb.codeocean.com/shiny-server-trusty/ ubuntu main' \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        pandoc=1.19.2.4~dfsg-1build4 \
        shiny-server=1.5.12.933 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("ggplot2", "3.4.0")' \
    && Rscript -e 'remotes::install_version("reshape2", "1.4.4")' \
    && Rscript -e 'remotes::install_version("shiny", "1.7.2")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "AnnotationHub", \
        "BSgenome.Hsapiens.UCSC.hg19", \
        "ConsensusClusterPlus", \
        "DESeq2", \
        "GenomicRanges", \
        "Repitools", \
        "TxDb.Hsapiens.UCSC.hg19.knownGene", \
        "annotatr", \
        "apeglm", \
        "clusterProfiler", \
        "enrichplot", \
        "fgsea", \
        "ggtree", \
        "org.Hs.eg.db" \
    ))' # Original versions: 3.6.0 1.4.3 1.62.0 1.38.1 1.50.1 1.44.0 3.2.2 1.24.0 1.20.0 4.6.0 1.18.1 1.24.0 3.6.2 3.16.0

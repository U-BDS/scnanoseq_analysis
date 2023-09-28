FROM bioconductor/bioconductor_docker:RELEASE_3_17

LABEL maintainer="Austyn Trull <atrull@uab.edu> and Lara Ianov <lianov@uab.edu>"
LABEL description="Environment which contains the tools used for evaluating the scnanoseq results"

RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libboost-all-dev \
    libxml2-dev \
    openjdk-8-jdk \
    python3-dev \
    python3-pip \
    wget \
    git \
    libfftw3-dev \
    libgsl-dev \
    llvm-11 \
    libbz2-dev \
    liblzma-dev \
    g++ \
    gcc

# CRAN Packages
RUN R --no-restore --no-save -e \
    "install.packages(c('remotes', \
                        'ggplot2', \
                        'dplyr', \
                        'grid', \
                        'sctransform', \
                        'clustree', \
                        'tidyverse', \
                        'SoupX', \
                        'rmarkdown', \
                        'markdown', \
                        'harmony', \
                        'hdf5r', \
                        'pheatmap', \
                        'ComplexUpset', \
                        'adabag', \
                        'gbm', \
                        'gprofiler2', \
                        'rsvd'))"

# Bioconductor pacakges
# as new datasets are released, discovered, some of the additional pacakges
# below will allows us to do more sophisticated analysis for validational purposes
RUN R --no-restore --no-save -e \
    "BiocManager::install(c('rtracklayer', \
                            'AnnotationHub', \
                            'biomaRt', \
                            'limma', \
                            'DESeq2', \
                            'DEXSeq', \
                            'EnhancedVolcano', \
                            'glmGamPoi', \
                            'multtest', \
                            'scater'))"                   

# Manually install Seurat, not wishing to use use v4.2.0
RUN R --no-restore --no-save -e "remotes::install_version('Seurat', version='4.3.0',repos='https://cloud.r-project.org', \
dependencies=c('Depends', 'Imports', 'LinkingTo'), upgrade='never')"

RUN R --no-restore --no-save -e "remotes::install_version('Matrix.utils', version='0.9.8',repos='https://cloud.r-project.org', \
dependencies=c('Depends', 'Imports', 'LinkingTo'), upgrade='never')"

# GitHub packages
RUN R --no-restore --no-save -e "devtools::install_github('chris-mcginnis-ucsf/DoubletFinder', upgrade='never')"

RUN R --no-restore --no-save -e "devtools::install_github('kostkalab/scds', ref='master', upgrade='never')"

RUN R --no-restore --no-save -e "devtools::install_github('13308204545/Chord', upgrade='never')"

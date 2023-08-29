FROM bioconductor/bioconductor_docker:RELEASE_3_17

LABEL maintainer="Austyn Trull <atrull@uab.edu>"
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

RUN mkdir -p /home/rstudio/scnanoseq_analysis/

# CRAN Packages
RUN R --no-restore --no-save -e \
    "install.packages(c('remotes', \
                        'ggplot2', \
                        'dplyr', \
                        'grid', \
                        'sctransform', \
                        'clustree'))"

# Manually install Seurat, not wishing to use use v4.2.0
RUN R --no-restore --no-save -e "remotes::install_version('Seurat', version='4.3.0',repos='https://cloud.r-project.org', \
dependencies=c('Depends', 'Imports', 'LinkingTo'), upgrade='never')"

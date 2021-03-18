FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base r-cran-randomforest python3.6 python3-pip python3-setuptools python3-dev

WORKDIR /app
RUN apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev

COPY requirements.txt /app/requirements.txt

RUN pip3 install -r requirements.txt

RUN Rscript -e "install.packages('data.table');install.packages('BiocManager', repos = 'https://cloud.r-project.org');BiocManager::install(c( 'devtools' , 'deconstructSigs' , 'MutationalPatterns' , 'remotes', 'data.table', 'dplyr', 'purrr', 'tidyr', 'furrr', 'Rcpp', 'cowplot', 'NMF', 'ggpubr', 'cli', 'reticulate', 'roxygen2'))"
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "library(devtools)"

RUN Rscript -e "library(reticulate);use_python('/usr/bin/python3.8');library(devtools);install_github('AlexandrovLab/SigProfilerMatrixGeneratorR');library(SigProfilerMatrixGeneratorR);install('GRCh37', rsync=FALSE, bash=TRUE)"
RUN R -e "BiocManager::install('BSgenome')" && \
    R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')" && \
    R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"
RUN R -e "BiocManager::install('ShixiangWang/sigminer@v1.2.3', dependencies = TRUE)"

COPY . /app



# FROM r-base:4.0.2


# LABEL \
#     author="Palash Pandey" \
#     maintainer="Palash Pandey" \
#     email="pp535@drexel.edu" \
#     description="Docker Image for MetaMutationalSigs" \
#     org.label-schema.license="Academic Free License v.3.0" \
#     org.label-schema.vcs-url="https://github.com/PalashPandey/MetaMutationalSigs"

# ## Install system dependencies
# RUN apt update -y && apt install -y libcurl4-openssl-dev libxml2-dev libssl-dev &&  \
#     apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*
# ## Install R packages which are easy to install
# RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')" && \
#     R -e "BiocManager::install(c('deconstructSigs' , 'MutationalPatterns' , 'remotes', 'data.table', 'dplyr', 'purrr', 'tidyr', 'furrr', 'Rcpp', 'cowplot', 'NMF', 'ggpubr', 'cli', 'reticulate', 'roxygen2'))"
# ## Install reference genome packages which are big
# RUN R -e "BiocManager::install('BSgenome')" && \
#     R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')" && \
#     R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"
# ## Install sigminer & SigProfiler & test SigProfiler
# RUN R -e "BiocManager::install('ShixiangWang/sigminer@v1.2.3', dependencies = TRUE)" && \
#     rm -rf /tmp/* /var/tmp/* && \
#     R -e "library(MutationalPatterns);library(deconstructSigs);library(sigfit);library('sigminer');data("cosmic_signatures_v2");data("cosmic_signatures_v3");load("data/signatures.exome.cosmic.v3.may2019.rda");load("data/signatures.cosmic.rda");load(system.file('extdata', 'toy_copynumber_tally_M.RData', package = 'sigminer', mustWork = TRUE)); mat = cn_tally_M[['nmf_matrix']]; print(mat); sigprofiler_extract(mat, '/opt/test_sp_install', range = 3:4, nrun = 2L, use_conda = TRUE); cat(paste(list.files('/opt/test_sp_install', recursive = TRUE), '\n')); sigprofiler_import('/opt/test_sp_install')" && \
#     rm -rf /opt/test_sp_install && \
#     /root/.local/share/r-miniconda/bin/conda clean -afy
# RUN pip install -r requirements.txt
# COPY . /usr/local/src/myscripts
# RUN R --vanilla -f /opt/pkg_check.R && \
#     R -e "install.packages('docopt', lib = .libPaths()[2])" && \
#     chmod u+x /opt/sigflow.R && ln -s /opt/sigflow.R /usr/bin/sigflow && \
#     cd /opt/test && chmod u+x test.sh && ./test.sh && rm -rf test_results && cd /root
# WORKDIR /usr/local/src/myscripts
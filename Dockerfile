FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base r-cran-randomforest python3.6 python3-pip python3-setuptools python3-dev


RUN apt-get -y install libcurl4-openssl-dev libssl-dev libxml2-dev

RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*


ADD ./ /app 
WORKDIR /app

RUN pip3 install -r requirements.txt
RUN python3.8 install_sigprofilerMatrixGenerator.py

RUN Rscript -e "install.packages('data.table');install.packages('BiocManager', repos = 'https://cloud.r-project.org');BiocManager::install(c( 'devtools' , 'deconstructSigs' , 'MutationalPatterns' , 'remotes', 'data.table', 'dplyr', 'purrr', 'tidyr', 'furrr', 'Rcpp', 'cowplot', 'NMF', 'ggpubr', 'cli', 'reticulate', 'roxygen2'))"
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "library(devtools)"

# RUN R -e "BiocManager::install('BSgenome')" && \
#     R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')" 

RUN R -e "install.packages('sigminer', dependencies = TRUE)"

RUN Rscript -e "library(reticulate);use_python('/usr/bin/python3.8');library(devtools);install_github('AlexandrovLab/SigProfilerMatrixGeneratorR');library(SigProfilerMatrixGeneratorR)"

# ENTRYPOINT ["python3.8" , "metaMutatationalSignatures.py"]
# CMD [ "--help" ]
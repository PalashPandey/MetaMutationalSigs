FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

ENV DOWNLOAD_STATIC_LIBV8 1

# RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3.6 python3-pip python3-setuptools python3-dev && \
#     apt update && apt-get -y install software-properties-common libcurl4-openssl-dev libssl-dev libxml2-dev libnode-dev vim-tiny wget && \

# apt-get purge -y r-base* r-recommended r-cran-*
# apt autoremove -y 
# apt update -y 


# DOWNLOAD_STATIC_LIBV8=1
ADD ./ /app 
WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3.6 python3-pip python3-setuptools python3-dev && \
    apt update && apt-get -y install software-properties-common libcurl4-openssl-dev libssl-dev libxml2-dev libnode-dev vim-tiny wget && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    apt update && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' &&\
    gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    apt install -y r-base r-base-core r-recommended r-base-dev && \
    add-apt-repository -y ppa:c2d4u.team/c2d4u4.0+ && apt-get install -f -y r-cran-ggplot2 r-cran-devtools r-cran-roxygen2 && \
    Rscript -e "update.packages(ask = FALSE, checkBuilt = TRUE);devtools::install_github('kgori/sigfit',build_opts = c('--no-resave-data', '--no-manual'));install.packages('BiocManager');BiocManager::install(dependencies = TRUE , c('MutationalPatterns', 'deconstructSigs' , 'ShixiangWang/sigminer@v2.0.0' ,  'dplyr', 'tidyr', 'ggpubr'))" && \
    pip3 install -r requirements.txt && python3.8 install_sigprofilerMatrixGenerator.py



# update.packages(ask = FALSE, checkBuilt = TRUE)
# install.packages("V8")
# install.packages("devtools", dependencies= TRUE)

# RUN Rscript -e "update.packages(ask = FALSE, checkBuilt = TRUE);install.packages("V8");devtools::install_github("kgori/sigfit",build_opts = c("--no-resave-data", "--no-manual"));install.packages('BiocManager');BiocManager::install(dependencies = TRUE , c('MutationalPatterns', 'deconstructSigs' , 'ShixiangWang/sigminer@v2.0.0' ,  'dplyr', 'tidyr', 'ggpubr'))"

EXPOSE 5001
ENTRYPOINT ["python3.8" , "metaMutatationalSignatures.py"]
# CMD [ "metaMutatationalSignatures.py" ]
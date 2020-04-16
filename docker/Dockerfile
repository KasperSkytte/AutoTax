FROM ubuntu:18.04

WORKDIR /opt

#install system dependencies using APT
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
  apt-get -y install --no-install-recommends --no-install-suggests \
    ca-certificates software-properties-common gnupg2 gnupg1 git \
    libcurl4-openssl-dev libssl-dev libxml2-dev wget parallel=20161222-1 && \
  mkdir -p ~/.parallel/will-cite #citing is out of context here, messes up TAP output

#add R CRAN mirror to apt sources to install R using APT and install R packages in parallel
#specific version of R to be installed
ENV R_BASE_VERSION 3.6.0
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
  add-apt-repository 'deb http://mirrors.dotsrc.org/cran/bin/linux/ubuntu bionic-cran35/' && \
  apt-get update && \
  apt-get -y install --no-install-recommends --no-install-suggests \
    r-base=${R_BASE_VERSION}* \
    r-base-dev=${R_BASE_VERSION}* \
    r-recommended=${R_BASE_VERSION}* && \
  mkdir -p ~/.R && \
  echo "MAKEFLAGS = -j100" > ~/.R/Makevars && \
  R -e "install.packages('BiocManager')" && \
  R -e "BiocManager::install(version = '3.9', ask = FALSE, Ncpus = 100)" && \
  R -e "BiocManager::install(c('Biostrings', 'doParallel', 'stringr', 'data.table', 'tidyr', 'dplyr'), Ncpus = 100, version = '3.9', ask = FALSE)"

#install SINA
RUN wget -q https://github.com/epruesse/SINA/releases/download/v1.6.0/sina-1.6.0-linux.tar.gz && \
  tar -zxf sina-1.6.0-linux.tar.gz && \
  rm sina-1.6.0-linux.tar.gz

#install BATS for unit testing
RUN git clone https://github.com/bats-core/bats-core.git && \
  /opt/bats-core/install.sh /usr/local

#make sure everything is in PATH
WORKDIR /autotax
ENV PATH="/autotax:/opt/sina-1.6.0-linux/bin:${PATH}"

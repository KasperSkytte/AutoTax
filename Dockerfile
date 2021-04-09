FROM ubuntu:18.04

WORKDIR /opt

### install system dependencies using APT
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
  apt-get -y install --no-install-recommends --no-install-suggests \
    ca-certificates \
    software-properties-common \
    gnupg2 \
    gnupg1 \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    wget \
    locales \
    parallel=20161222-1 && \
  mkdir -p ~/.parallel/will-cite #citing is out of context here, messes up TAP output

### generate and set up locales
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
  locale-gen en_US.utf8 && \
  /usr/sbin/update-locale LANG=en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

### install and setup R and the R library
#specific version of R to be installed
ENV R_BASE_VERSION 3.6.0

#add R CRAN mirror to apt sources to install R using APT,
#install R packages in parallel, use specific AAU CRAN mirror,
#remove user library from .libPaths() as it may be used instead
#if user has installed any pkgs on host system and $HOME is mounted
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
  add-apt-repository 'deb http://mirrors.dotsrc.org/cran/bin/linux/ubuntu bionic-cran35/' && \
  apt-get update && \
  apt-get -y install --no-install-recommends --no-install-suggests \
    r-base=${R_BASE_VERSION}* \
    r-base-dev=${R_BASE_VERSION}* \
    r-recommended=${R_BASE_VERSION}* && \
  mkdir -p ~/.R && \
  echo "MAKEFLAGS = -j100" > ~/.R/Makevars && \
  echo "options(repos = c(CRAN = 'https://mirrors.dotsrc.org/cran'), download.file.method = 'libcurl')" >> /usr/lib/R/etc/Rprofile.site && \
  sed -i s/^R_LIBS_USER=/#R_LIBS_USER=/g /etc/R/Renviron &&\
  R -e "install.packages('BiocManager')" && \
  R -e "BiocManager::install(version = '3.9', ask = FALSE, Ncpus = 100)" && \
  R -e "BiocManager::install(c('Biostrings', 'doParallel', 'stringr', 'data.table', 'tidyr', 'dplyr'), Ncpus = 100, version = '3.9', ask = FALSE)"

### install SINA
RUN wget -q https://github.com/epruesse/SINA/releases/download/v1.6.0/sina-1.6.0-linux.tar.gz && \
  tar -zxf sina-1.6.0-linux.tar.gz && \
  rm sina-1.6.0-linux.tar.gz


### install vsearch
RUN wget -q https://github.com/torognes/vsearch/releases/download/v2.17.0/vsearch-2.17.0-linux-x86_64.tar.gz && \
  tar -zxf vsearch-2.17.0-linux-x86_64.tar.gz && \
  rm vsearch-2.17.0-linux-x86_64.tar.gz

### install BATS for unit testing
RUN git clone https://github.com/bats-core/bats-core.git && \
  /opt/bats-core/install.sh /usr/local

### copy AutoTax repo into /opt/autotax
COPY . /opt/autotax/
RUN chmod +x /opt/autotax/autotax.bash

### make sure everything is in PATH
ENV PATH="/autotax:/opt/sina-1.6.0-linux/bin/:/opt/vsearch-2.17.0-linux-x86_64/bin:${PATH}"

WORKDIR /autotax
ENTRYPOINT ["bash", "/opt/autotax/autotax.bash"]
CMD ["-h"]

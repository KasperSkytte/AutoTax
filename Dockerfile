FROM ubuntu:bionic-20210930

WORKDIR /opt

### install system dependencies using APT
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update -qqy \
  && apt-get -y install --no-install-recommends --no-install-suggests \
    ca-certificates \
    software-properties-common \
    gnupg2 \
    gnupg1 \
    git \
    wget \
    locales \
    parallel=20161222-1 \
  && mkdir -p ~/.parallel/will-cite #citing is out of context here, messes up TAP output

### generate and set up locales
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen en_US.utf8 \
  && /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

### install and setup R and the R library
#specific version of R to be installed
ENV R_BASE_VERSION 3.6.3

#install R, remove user library from .libPaths() as it may be used instead
#if user has installed any pkgs on host system and $HOME is mounted
RUN apt-get -qqy update \
  && DEBIAN_FRONTEND=noninteractive apt-get -y install --no-install-recommends --no-install-suggests \
    gdebi-core \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    libcairo2-dev \
    pandoc \
    curl \
  && curl -O https://cdn.rstudio.com/r/ubuntu-1804/pkgs/r-${R_BASE_VERSION}_1_amd64.deb \
  && gdebi --non-interactive r-${R_BASE_VERSION}_1_amd64.deb \
  && rm r-${R_BASE_VERSION}_1_amd64.deb \
  && ln -s /opt/R/${R_BASE_VERSION}/bin/R /usr/local/bin/R \
  && ln -s /opt/R/${R_BASE_VERSION}/bin/Rscript /usr/local/bin/Rscript \
  && mkdir -p ~/.R \
  && echo "MAKEFLAGS = -j" > ~/.R/Makevars \
  && echo "options(repos = c(CRAN = 'https://packagemanager.rstudio.com/all/__linux__/bionic/2021-11-17+MTo2MzM5NDY2LDI6NDUyNjIxNTs1MDMxRDg3Qg'), download.file.method = 'libcurl')" >> /opt/R/${R_BASE_VERSION}/lib/R/etc/Rprofile.site \
  && sed -i s/^R_LIBS_USER=/#R_LIBS_USER=/g /opt/R/${R_BASE_VERSION}/lib/R/etc/Renviron \
  && R -e "install.packages('BiocManager')" \
  && R -e "BiocManager::install(version = '3.9', ask = FALSE, Ncpus = 100)" \
  && R -e "BiocManager::install(c('Biostrings', 'doParallel', 'stringr', 'data.table', 'tidyr', 'dplyr'), Ncpus = 100, version = '3.9', ask = FALSE)"

### install SINA
RUN wget -q https://github.com/epruesse/SINA/releases/download/v1.6.0/sina-1.6.0-linux.tar.gz \
  && tar -zxf sina-1.6.0-linux.tar.gz \
  && rm sina-1.6.0-linux.tar.gz

### install usearch 32-bit limited free version
RUN wget -q http://drive5.com/downloads/usearch11.0.667_i86linux32.gz -O usearch11.gz \
  && gunzip usearch11.gz \
  && chmod +x usearch11 \
  && mkdir -p usearch11.0.667_i86linux32 \
  && mv usearch11 -t usearch11.0.667_i86linux32/
  

### install vsearch
RUN wget -q https://github.com/torognes/vsearch/releases/download/v2.17.0/vsearch-2.17.0-linux-x86_64.tar.gz \
  && tar -zxf vsearch-2.17.0-linux-x86_64.tar.gz \
  && rm vsearch-2.17.0-linux-x86_64.tar.gz

### install BATS for unit testing
RUN wget -q https://github.com/bats-core/bats-core/archive/refs/tags/v1.3.0.tar.gz \
  && tar -zxf v1.3.0.tar.gz \
  && rm v1.3.0.tar.gz \
  && /opt/bats-core-1.3.0/install.sh /usr/local \
  && rm -rf /opt/bats-core-1.3.0

### copy AutoTax repo into /opt/autotax
COPY . /opt/autotax/
RUN chmod +x /opt/autotax/autotax.bash

### make sure everything is in PATH
ENV PATH="/opt/usearch11.0.667_i86linux32:${PATH}"
ENV PATH="/opt/sina-1.6.0-linux/bin:${PATH}"
ENV PATH="/opt/vsearch-2.17.0-linux-x86_64/bin:${PATH}"
ENV PATH="/opt/autotax:${PATH}"
ENV PATH="/autotax:${PATH}"

### clean up
RUN rm -rf /var/lib/apt/lists/* /var/cache/apt/*

VOLUME /autotax
WORKDIR /autotax
ENTRYPOINT ["bash", "/opt/autotax/autotax.bash"]
CMD ["-h"]

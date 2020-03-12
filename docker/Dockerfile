FROM ubuntu:18.04

WORKDIR /opt

#install system dependencies using APT
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
RUN apt-get -y install --no-install-recommends --no-install-suggests \
        ca-certificates software-properties-common gnupg2 gnupg1 git \
        parallel=20161222-1 libcurl4-openssl-dev libssl-dev libxml2-dev wget

#add R CRAN mirror to apt sources to install R using APT
#specific version of R to be installed
ENV R_BASE_VERSION 3.6.0
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 
RUN add-apt-repository 'deb http://mirrors.dotsrc.org/cran/bin/linux/ubuntu bionic-cran35/'
RUN apt-get update
RUN apt-get -y install r-base=${R_BASE_VERSION}* \
		r-base-dev=${R_BASE_VERSION}* \
		r-recommended=${R_BASE_VERSION}*

#R packages, use latest version from Bioconductor release 3.9
#use all processors available with Ncpus = 100 to install faster
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.9', ask = FALSE, Ncpus = 100)"
RUN R -e "BiocManager::install(c('Biostrings', 'doParallel', 'stringr', 'data.table', 'tidyr', 'dplyr'), Ncpus = 100, version = '3.9', ask = FALSE)"

#install SINA
RUN wget -q https://github.com/epruesse/SINA/releases/download/v1.6.0/sina-1.6.0-linux.tar.gz
RUN tar -zxf sina-1.6.0-linux.tar.gz
RUN rm sina-1.6.0-linux.tar.gz

#install BATS for unit testing
RUN git clone https://github.com/bats-core/bats-core.git
RUN /opt/bats-core/install.sh /usr/local
#silence citation notice from GNU parallel (not used in a context where it makes sense to cite)
RUN mkdir -p ~/.parallel && touch ~/.parallel/will-cite

#make sure everything is in PATH
WORKDIR /home/user/
ENV PATH="/home/user:/opt/sina-1.6.0-linux/bin:${PATH}"

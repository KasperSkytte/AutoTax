FROM ubuntu:18.04

WORKDIR /opt
ENV DEBIAN_FRONTEND noninteractive
ENV R_BASE_VERSION 3.6.0

RUN apt-get update
RUN apt-get install -y wget
#install R
RUN apt-get -y install --no-install-recommends --no-install-suggests \
        ca-certificates software-properties-common gnupg2 gnupg1 
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 
RUN add-apt-repository 'deb http://mirrors.dotsrc.org/cran/bin/linux/ubuntu bionic-cran35/'
RUN apt-get update
RUN apt-get -y install r-base=${R_BASE_VERSION}* \
		r-base-dev=${R_BASE_VERSION}* \
		r-recommended=${R_BASE_VERSION}*

#install SINA
RUN wget -q https://github.com/epruesse/SINA/releases/download/v1.6.0/sina-1.6.0-linux.tar.gz
RUN tar -zxf sina-1.6.0-linux.tar.gz
RUN rm sina-1.6.0-linux.tar.gz

#make sure everything is in PATH
RUN echo "export PATH=$PATH:/opt/sina-1.6.0-linux/bin" >> /root/.bashrc
WORKDIR /home/user/

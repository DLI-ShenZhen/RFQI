#################################################
#### Result image is :v1.1 ########
#### Based on ubuntu:16.04 #####################
#################################################

FROM ubuntu:16.04
MAINTAINER luanenhui luanenhui2009@gmail.com

RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list
RUN apt-key adv --recv-keys --keyserver keyserver.ubuntu.com 51716619E084DAB9  # secure APT source for r-project
RUN apt update
RUN apt install -y apt-transport-https

RUN apt upgrade -y

RUN apt install -y gcc cmake libcurl14-openssl-dev libxml2-dev curl libcairo2-dev libssl-dev vim wget git
RUN apt install -y r-base

RUN git clone https://github.com/luanenhui/RFQI.git
WRRKDIR RFQI

########### install dependency ###########################
RUN tar -xzf dependency/hdf5-1.8.18.tar.gz
RUN tar -xzf dependency/zlib-1.2.8.tar.gz
RUN tar -xzf dependency/netcdf-4.4.1.1.tar.gz

RUN cd zlib-1.2.8 && ./configure --prefix=/usr/local && make install && cd ../ && rm -rf zlib-1.2.8
RUN cd hdf5-1.8.18 && ./configure --with-zlib=/usr/local --prefix=/usr/local && make install && cd ../ && rm -rf hdf5-1.8.18
RUN cd netcdf-4.4.1.1 && CPPFLAGS=-I/usr/local/include LDFLAGS=-L/usr/local/lib ./configure --prefix=/usr/local && make install && cd ../ && rm -rf netcdf-4.4.1.1

########### install dependency R packages ###########################
RUN Rscript -e 'install.packages(c("devtools", "future", "scales"), repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN")'
RUN Rscript -e 'options(download.file.method="wget");source("https://bioconductor.org/biocLite.R"); biocLite(suppressUpdates=TRUE)'
RUN Rscript -e 'options(download.file.method="wget");library(BiocInstaller); biocLite("mzR",suppressUpdates=TRUE)'
RUN Rscript -e 'options(download.file.method="wget");library(BiocInstaller); biocLite("MSnbase",suppressUpdates=TRUE)'
RUN Rscript -e 'options(download.file.method="wget");library(BiocInstaller); biocLite("xcms",suppressUpdates=TRUE)'

########### install dependency R packages ###########################
RUN Rscript -e 'options(download.file.method="wget";library(devtools);devtools::install("RFQI")'

###### delete RFQI and change workdir
RUN rm -rf RFQI
WORKDIR /

ENTRYPOINT ["/bin/bash"]


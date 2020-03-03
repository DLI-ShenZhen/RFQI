#################################################
#### Result image is lipid_analyzer:v1.0 ########
#### Based on ubuntu:latest #####################
#################################################

FROM ubuntu:16.04
MAINTAINER luanenhui luanenhui@icarbonx.com

RUN apt update
RUN apt install -y apt-transport-https
ADD sources.list /etc/apt/

RUN apt-key adv --recv-keys --keyserver keyserver.ubuntu.com 51716619E084DAB9
RUN apt-key adv --recv-keys --keyserver keyserver.ubuntu.com D259B7555E1D3C58
RUN apt-key adv --recv-keys --keyserver keyserver.ubuntu.com 7EA0A9C3F273FCD8
RUN apt update
RUN apt upgrade -y

RUN apt install -y gcc cmake r-base curl libcairo2-dev libssl-dev vim wget
RUN mkdir /home/software

########### install dependency ###########################
# ADD will auto untar tar.gz files
Add dependency/zlib-1.2.8.tar.gz /usr/local/bin/
Add dependency/hdf5-1.8.18.tar.gz /usr/local/bin/
Add dependency/netcdf-4.4.1.1.tar.gz /usr/local/bin/

RUN cd /usr/local/bin/zlib-1.2.8 && ./configure --prefix=/usr/local && make install
RUN cd /usr/local/bin/hdf5-1.8.18 && ./configure --with-zlib=/usr/local --prefix=/usr/local && make install
RUN cd /usr/local/bin/netcdf-4.4.1.1 && CPPFLAGS=-I/usr/local/include LDFLAGS=-L/usr/local/lib ./configure --prefix=/usr/local && make install

######### set LD_LIBRARY_PATH #######################
ENV LD_LIBRARY_PATH /usr/local/lib

######### install gocos ##################################
Add dependency/gocos /usr/local/bin
RUN chmod u+x /usr/local/bin/gocos
Add dependency/cos.config.json /root/.cos.config.json

######## install R packages #############################

COPY dependency/MetMatch_1.0.2.tar.gz /home/software/MetMatch_1.0.2.tar.gz
COPY dependency/DDADataProcesser_1.0.4.tar.gz /home/software/DDADataProcesser_1.0.4.tar.gz
COPY dependency/LipidAnalyzer_3.0.1.tar.gz /home/software/LipidAnalyzer_3.0.1.tar.gz

########## install R packages ##########################
RUN apt install -y libxml2-dev
RUN Rscript -e 'install.packages("Rcpp", repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN")'
RUN Rscript -e 'options(download.file.method="wget");source("https://bioconductor.org/biocLite.R"); biocLite(suppressUpdates=TRUE)'
RUN Rscript -e 'options(download.file.method="wget");library(BiocInstaller); biocLite("BiocParallel",suppressUpdates=TRUE)'
RUN Rscript -e 'options(download.file.method="wget");library(BiocInstaller); biocLite("igraph",suppressUpdates=TRUE)'
#RUN Rscript -e 'library(BiocInstaller); biocLite("mzR", suppressUpdates=TRUE)'
RUN Rscript -e 'options(download.file.method="wget");library(BiocInstaller); biocLite("xcms",suppressUpdates=TRUE)'
RUN Rscript -e 'options(download.file.method="wget");library(BiocInstaller); biocLite("CAMERA",suppressUpdates=TRUE)'
RUN Rscript -e 'install.packages("/home/software/MetMatch_1.0.2.tar.gz", repos=NULL, type="source")'
RUN Rscript -e 'install.packages("/home/software/DDADataProcesser_1.0.4.tar.gz", repos=NULL, type="source")'

COPY dependency/MetAnalyzer_3.2.0.tar.gz /home/software/MetAnalyzer_3.2.0.tar.gz
RUN Rscript -e 'install.packages("/home/software/MetAnalyzer_3.2.0.tar.gz", repos=NULL, type="source")'


#RUN apt install -y apt-file
#RUN apt-file update
#RUN apt-file find Intrinsic.h  # find libxt-dev including Intrinsic.h
RUN apt install -y libxt-dev
RUN Rscript -e 'install.packages("Cairo", repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN")'

Add src/run_MetCleaning.pl /usr/local/bin/
Add src/run_MetCleaning.r /usr/local/bin/
Add shell/met_analyze.sh /usr/local/bin/met_analyze.sh
RUN chmod u+x /usr/local/bin/met_analyze.sh

ENTRYPOINT ["/usr/local/bin/met_analyze.sh"]


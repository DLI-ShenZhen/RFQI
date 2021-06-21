#################################################
#### Result image is :v1.1 ########
#### Based on ubuntu:18.04 #####################
#################################################

FROM ubuntu:18.04
MAINTAINER luanenhui luanenhui2009@gmail.com

RUN chmod 777 /tmp
RUN apt update
RUN apt install -y gnupg gcc cmake libcurl4-openssl-dev libxml2-dev curl libcairo2-dev libssl-dev vim wget git
RUN export DEBIAN_FRONTEND=noninteractive && apt install -y tzdata
RUN apt upgrade -y

# revise timezone to Asia/Shanghai
ENV TZ=Asia/Shanghai

# install r-base
RUN apt install -y --no-install-recommends software-properties-common dirmngr
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN apt install -y --no-install-recommends r-base r-base-dev

# install RFQI
RUN git clone https://github.com/luanenhui/RFQI.git
WORKDIR RFQI

########### install dependency ###########################
RUN tar -xzf dependency/hdf5-1.8.18.tar.gz
RUN tar -xzf dependency/zlib-1.2.8.tar.gz
RUN tar -xzf dependency/netcdf-4.4.1.1.tar.gz

RUN cd zlib-1.2.8 && ./configure --prefix=/usr/local && make install && cd ../ && rm -rf zlib-1.2.8
RUN cd hdf5-1.8.18 && ./configure --with-zlib=/usr/local --prefix=/usr/local && make install && cd ../ && rm -rf hdf5-1.8.18
RUN apt install -y m4 && cd netcdf-4.4.1.1 && CPPFLAGS=-I/usr/local/include LDFLAGS=-L/usr/local/lib ./configure --prefix=/usr/local && make install && cd ../ && rm -rf netcdf-4.4.1.1

########### install dependency R packages ###########################
RUN Rscript -e 'install.packages("BiocManager", repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN")'
RUN Rscript -e 'BiocManager::install(version = "3.13", update=FALSE)'
RUN Rscript -e 'BiocManager::install(c("mzR", "MSnbase", "igraph", "xcms"), update=FALSE)'
RUN Rscript -e 'BiocManager::install(c("devtools"), update=FALSE)'

########### install dependency R packages ###########################
WORKDIR /
RUN Rscript -e 'library(devtools);devtools::install("/RFQI")'

###### delete RFQI and change workdir
RUN rm -rf RFQI

ENTRYPOINT ["/bin/bash"]


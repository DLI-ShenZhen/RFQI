### Install dependency 
> RFQI depends on R-3.4, R-3.6 is not suitable 

```bash
# install R 3.4.4
sudo echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list
sudo apt-key adv --recv-keys --keyserver keyserver.ubuntu.com 51716619E084DAB9  # secure APT source for r-project
sudo apt update
sudo apt install libcurl4-openssl-dev curl libssl-dev r-base

# install netcdf4
tar -xzf dependency/hdf5-1.8.18.tar.gz
tar -xzf dependency/zlib-1.2.8.tar.gz
tar -xzf dependency/netcdf-4.4.1.1.tar.gz

# The flowing codes need sudoer priviledge 
cd zlib-1.2.8 && ./configure --prefix=/usr/local && make install && cd ../ && rm -rf zlib-1.2.8
cd hdf5-1.8.18 && ./configure --with-zlib=/usr/local --prefix=/usr/local && make install && cd ../ && rm -rf hdf5-1.8.18
cd netcdf-4.4.1.1 && CPPFLAGS=-I/usr/local/include LDFLAGS=-L/usr/local/lib ./configure --prefix=/usr/local && make install && cd ../ && rm -rf netcdf-4.4.1.1

# install rstudio server
sudo apt-get install gdebi-core
wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-1.2.5033-amd64.deb
sudo gdebi rstudio-server-1.2.5033-amd64.deb
```

### Install dependency R packages
```r
options(repos = "https://mirrors.ustc.edu.cn/CRAN/")
install.package("devtools")

source("https://bioconductor.org/biocLite.R")
biocLite("mzR=2.16.1")
biocLite("MSnbase=2.8.3")
biocLite("xcms=3.4.4")

devtools::install_github("https://github.com/luanenhui/RFQI.git")
```

### [tutorial](./doc/tutorial.html)

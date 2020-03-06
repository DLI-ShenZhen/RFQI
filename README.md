RFQI is depends on R-3.4.4 OR R-3.5.3;  R-3.6 is not suitable
---
## Installation
### Install R-3.4.4 and dependency packages
#### For linux 
* ubuntu 16.04 can install R-3.4.4

```bash
# install R 3.4.4
sudo echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list
sudo apt-key adv --recv-keys --keyserver keyserver.ubuntu.com 51716619E084DAB9  # secure APT source for r-project
sudo apt update
sudo apt install -y libcurl4-openssl-dev curl libssl-dev libxml2-dev
sudo apt install -y r-base 

# install netcdf4
tar -xzf dependency/hdf5-1.8.18.tar.gz
tar -xzf dependency/zlib-1.2.8.tar.gz
tar -xzf dependency/netcdf-4.4.1.1.tar.gz

# The flowing codes need sudoer priviledge 
cd zlib-1.2.8 && ./configure --prefix=/usr/local && make install && cd ../ && rm -rf zlib-1.2.8
cd hdf5-1.8.18 && ./configure --with-zlib=/usr/local --prefix=/usr/local && make install && cd ../ && rm -rf hdf5-1.8.18
cd netcdf-4.4.1.1 && CPPFLAGS=-I/usr/local/include LDFLAGS=-L/usr/local/lib ./configure --prefix=/usr/local && make install && cd ../ && rm -rf netcdf-4.4.1.1
```

#### For windows
You need install both of r-base and Rtools  
[r-base 3.4.4](https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/windows/base/old/3.4.4/R-3.4.4-win.exe)  
[Rtools](https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/windows/Rtools/Rtools35.exe)  

---

### Install dependency R packages
```r
options(repos = "https://mirrors.ustc.edu.cn/CRAN/")
install.packages("devtools")
install.packages("future")

source("https://bioconductor.org/biocLite.R")
biocLite("mzR")
biocLite("MSnbase")
biocLite("xcms")

```
---
### Install RFQI from github and test
```r
devtools::install_github("https://github.com/luanenhui/RFQI.git")

```
---

### [vignette tutorial](./doc/tutorial.html)

### Release history
| version  | date | author  | Note  |
|-------|:---:|-----------|-------:|
| *v1.1*| 2020-03-05 | luanenhui2009@gmail.com | Replace mclapply with parLapply for windows |
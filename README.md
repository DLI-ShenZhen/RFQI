RFQI is depends on R-4.1.0 OR greater
---
## Installation
### Install R-4.1.0 and dependency packages
#### [For linux](https://mirrors.tuna.tsinghua.edu.cn/CRAN/) 

```bash
# install R 4.1.0
apt update -qq
apt install -y --no-install-recommends software-properties-common dirmngr
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
apt install -y --no-install-recommends r-base r-base-dev


# install netcdf4
tar -xzf dependency/hdf5-1.8.18.tar.gz
tar -xzf dependency/zlib-1.2.8.tar.gz
tar -xzf dependency/netcdf-4.4.1.1.tar.gz

# The flowing codes need sudoer priviledge 
cd zlib-1.2.8 && ./configure --prefix=/usr/local && make install && cd ../ && rm -rf zlib-1.2.8
cd hdf5-1.8.18 && ./configure --with-zlib=/usr/local --prefix=/usr/local && make install && cd ../ && rm -rf hdf5-1.8.18
cd apt install -y m4 && netcdf-4.4.1.1 && CPPFLAGS=-I/usr/local/include LDFLAGS=-L/usr/local/lib ./configure --prefix=/usr/local && make install && cd ../ && rm -rf netcdf-4.4.1.1
```

#### For windows
You need install both of r-base and Rtools  
[r-base 4.1.0](https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/windows/base/R-4.1.0-win.exe)  
[Rtools](https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/windows/Rtools/rtools40v2-x86_64.exe)  

---

### Install dependency R packages
```r
install.packages("BiocManager", repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN")
BiocManager::install(version = "3.13", update=FALSE)
BiocManager::install(c("mzR", "MSnbase", "igraph", "xcms"), update=FALSE)
BiocManager::install(c("devtools"), update=FALSE)

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
| *v2.0.0* | 2021-06-21 | luanenhui2009@gmail.com | switch from R 3.4.4 to R 4.1.0 |
| *v1.1*| 2020-03-05 | luanenhui2009@gmail.com | Replace mclapply with parLapply for windows |
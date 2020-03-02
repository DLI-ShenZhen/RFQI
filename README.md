### The first version of this package is called *XYMeta*, but later it will be changed into *RFQI*  

XYMeta depends on xcms

For R >= 3.6
```r
if (!requireNamespace("BiocManager", quitely=TRUE))
  intall.packages("BiocManager")

BiocManager::install("MSnbase")
BiocManager::install("mzR")
BiocManager::install("xcms")
```

For R < 3.6
```r
source("https://bioconductor.org/biocLite.R")
biocLite("MSnbase")
biocLite("mzR")
biocLite("xcms")
```

Then you can install xymeta

```r
# download xymeta_0.9.tar.gz to your computer
install.packages("xymeta_0.9.tar.gz", repos=NULL)

library(xymeta)
help(package="xymeta")
```

### [tutorial](./doc/tutorial.html)

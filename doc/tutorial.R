## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message=FALSE
)

## -----------------------------------------------------------------------------
library(RFQI)
options(stringsAsFactors = FALSE)

## ---- include=TRUE------------------------------------------------------------
files <- dir(system.file(package="RFQI", dir="extdata"), 
              full.name=TRUE, 
              pattern="mzXML$")
param = new("ParamSet", binSize=0.25, peakwidth=c(2,30), ppm=10, noise=0, absMz=0.005, absRt=15)

## -----------------------------------------------------------------------------
refFile = files[1]
files = files[2]

## -----------------------------------------------------------------------------
features = getPeakRange_multipleFiles(f.in = files, ref = refFile, param=param, cores=2)
# save(features, ref, files, file="test_features.rData")

## -----------------------------------------------------------------------------
str(features, max.level = 1)

## -----------------------------------------------------------------------------
# data("features", package="RFQI")

## -----------------------------------------------------------------------------
files <- dir(system.file(package="RFQI", dir="extdata"), full.name=TRUE, pattern="mzXML$")
MS2 = getMS2FromFiles(files=files, ref_file = refFile)
str(MS2, max.level = 1)

## -----------------------------------------------------------------------------
mzCount = unlist(lapply(MS2$MS2, nrow)) # A qualified MS2 spectrum should have at least 3 fragments
mzMaxInt = unlist(lapply(MS2$MS2, function(x){max(x[,2])}))  # A qualified MS2 spectrum should have at least one fragment whose intensity greater than 32 (user defined)
CountFlag = !(mzCount %in% c(0,1,2))
IntFlag = mzMaxInt > 50
flag = CountFlag & IntFlag

MS2$MS2 = MS2$MS2[flag]
MS2$precursorInfo = MS2$precursorInfo[flag,]

## -----------------------------------------------------------------------------
names(MS2$MS2) = rownames(MS2$precursorInfo) = sapply(1:nrow(MS2$precursorInfo),
convert_num_to_char, prefix = "M", n=nchar(nrow(MS2$precursorInfo)))

## -----------------------------------------------------------------------------
db.MS2 = align_MS2_to_MS1(ms2Info = MS2, features = features)

## -----------------------------------------------------------------------------
# we can calcalate similarity distribution for one feature
idx = getFeatureHasMS2(MS2DB=db.MS2, n=2)  # n is the minimum number of MS2 belonging to the same feature
distri = get_ms2Cor_inner(MS2DB = db.MS2, cores1=2, cores2=2, maxMS2 = 100, idx=idx[1])

## -----------------------------------------------------------------------------
plot(density(distri[[1]][,1], na.rm=TRUE), main="MS2 similarity distribution", col="red", xlab="cosine similarity score")
for (i in 2:nrow(distri[[1]])){
  lines(density(distri[[1]][,i], na.rm=TRUE), col="red")
}

## -----------------------------------------------------------------------------
db.MS2.subset = list(MS2=db.MS2$MS2[1:2000], MS2_to_MS1=db.MS2$MS2_to_MS1[1:2000,]) # this line is only for test
# db.MS2.subset = db.MS2
ms2InnerCor = get_ms2Cor_inner(MS2DB = db.MS2.subset,cores1 = 2, cores2=2, n=2, maxMS2 = 100)

## -----------------------------------------------------------------------------
data("spectrumDB", package = "RFQI")
str(spectrumDB, max.level = 1)

## -----------------------------------------------------------------------------
print("beginning identify")
anno = get_identify(lib = spectrumDB, MS2DB = db.MS2.subset, MS2_inner_cor = ms2InnerCor, absMz = 0.045, adduct="RP_pos", cores=3)
str(anno, max.level = 1)
print("end identify")

## -----------------------------------------------------------------------------
result = list()
for (file in files){
  result[[file]] = getPeaks_L(file = file, peak_range = features$group_FT, ref_file = refFile, MS1=TRUE, step=0.1, MS2=FALSE)$peaks[,"into",drop=FALSE]
}
intensity = do.call(cbind,result)
colnames(intensity) = names(result)
str(intensity, max.level = 1)


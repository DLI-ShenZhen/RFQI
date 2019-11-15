library(xymeta)
set.seed(22)
# First, we need to build a reference feature table using a set of mzXML files
mzXMLFiles = list.files(path="PaperData", pattern="mzXML$", full.names=TRUE, recursive = TRUE)
mzXMLFiles = mzXMLFiles[-grep(pattern="reference", mzXMLFiles)]

seedFiles = sample(mzXMLFiles, 85)

# we also need to set a reference mzXML
refFile = list.files(path="PaperData/MSS170428001_Urine/", pattern="reference", full.names=TRUE, recursive = TRUE)

# Then we extract peaks from mzXML files and group them
# First, we need to set parameter for retention time alignment(binSize), peak detection(ppm, noise, peakwidth) and peak group(absMz, absRt)
param = new("ParamSet", binSize=0.25, peakwidth=c(2,30), ppm=10, noise=0, absMz=0.005, absRt=15)

# Then we execute function *getPeakRange_multipleFiles*, this function extracts peaks from files, and group peaks to features
# This step is time-consuming, so we save the result for later use
features = getPeakRange_multipleFiles(f.in=seedFiles, ref=refFile, param=param, cores=2)
save(seedFiles, features, file="randomFileFeatureTable.rData")

# features contains original peaks and grouped features
str(features, max.level = 1)

# we can plot feature Area
plotArea(refPeaks = features$group_FT)

# ---- extract MS2 ---- #
# Next, we need to extract MS2 for annotating feature table
MS2 = getMS2FromFiles(sample(mzXMLFiles,10), ref_file = refFile)

# we can calculate some metrics for MS2 quality
mzCount = unlist(lapply(MS2$MS2, nrow))
mzMaxInt = unlist(lapply(MS2$MS2, function(x){max(x[,2])}))

flag = !(mzCount %in% c(0,1,2))
MS2$MS2 = MS2$MS2[flag]
MS2$precursorInfo = MS2$precursorInfo[flag,]

# we can set names for each MS2 
names(MS2$MS2) = rownames(MS2$precursorInfo) = sapply(1:nrow(MS2$precursorInfo), 
                                                      convert_num_to_char, prefix = "M", n=nchar(nrow(MS2$precursorInfo)))

# We need to map MS2 to feature based on precursor's m/z and rt
db.MS2 = align_MS2_to_MS1(ms2Info = MS2, features = features)
# save(db.MS2, file="randomDbMS2.rData")

# Next, we calculate MS2 similarity matrix belonging to same feature
# the time complexity is O(n2), where n represents the number of MS2 belonging the the same feature, we can limit the max number using maxMS2 parameter
head(sort(table(db.MS2$MS2_to_MS1[,"MS1"]), decreasing=TRUE))
ms2InnerCor = get_ms2Cor_inner(MS2DB = db.MS2,cores1 = 3, cores2 = 1, n=2, maxMS2 = 100)

# ---- identification ------ #


# ---- quantification ---- #
result = list()
for (file in mzXMLFiles){
  result[[file]] = getPeaks_L(file = file, peak_range = features$group_FT, ref_file = refFile, MS1=TRUE, step=0.1, MS2=FALSE)$peaks[,"into"]
}
intensity = do.call(cbind, result)
colnames(intensity) = sapply(basename(colnames(intensity)), function(x){strsplit(x, split=".", fixed = TRUE)[[1]][1]})
save(intensity, file="randomFTIntensity.rData")

# intensity = lapply(mzXMLFiles, getPeaks_L, peak_range = features$group_FT, ref_file = refFile, MS1=TRUE, step=0.1, MS2=FALSE)



# ----- check rFT quality by evaluating similarity matrix ---- #







# # and we can also plot peaks and features area to compare grouping process
# df = rbind(features$peaks[,c("mzmin", "mzmax","rtmin","rtmax")], features$group_FT[,c("mzmin","mzmax","rtmin","rtmax")])
# group = rep(c("peaks","features"), c(nrow(features$peaks), nrow(features$group_FT)))
# plotArea(refPeaks = df, group=group)
# 
# plotArea(refPeaks = df, group = group, mzr=c(400,500), rtr=c(180,220))


# select seed mzXML files to build
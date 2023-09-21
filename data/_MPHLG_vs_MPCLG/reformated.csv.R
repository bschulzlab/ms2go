library(MSstats)
ms<-read.csv("/data/_MPHLG_vs_MPCLG/reformated.csv", sep=",")
quant<-dataProcess(ms)
levels(quant$GROUP_ORIGINAL)
comparison0<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0),nrow=1)
comparison<-comparison0
row.names(comparison)<-c("(1)MPHLGvs(-1)MPCLG")
testResultMultiComparisons<-groupComparison(contrast.matrix=comparison,data=quant,labeled=FALSE,interference=FALSE,featureVar=TRUE)
write.csv(testResultMultiComparisons$ComparisonResult, "/data/_MPHLG_vs_MPCLG/reformated.csv_msstats.csv")

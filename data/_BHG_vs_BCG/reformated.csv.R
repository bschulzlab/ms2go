library(MSstats)
ms<-read.csv("/data/_BHG_vs_BCG/reformated.csv", sep=",")
quant<-dataProcess(ms)
levels(quant$GROUP_ORIGINAL)
comparison0<-matrix(c(-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
comparison<-comparison0
row.names(comparison)<-c("(1)BHGvs(-1)BCG")
testResultMultiComparisons<-groupComparison(contrast.matrix=comparison,data=quant,labeled=FALSE,interference=FALSE,featureVar=TRUE)
write.csv(testResultMultiComparisons$ComparisonResult, "/data/_BHG_vs_BCG/reformated.csv_msstats.csv")

library(GOstats)
library(GSEABase)
GODB<-read.table("/data/_LHG_vs_LCG/(1)LHGvs(-1)LCGgostats_association.txt", sep="	", header=TRUE, fill=TRUE, comment.char = "")

        universe<-read.table("/data/_LHG_vs_LCG/_universe.txt", sep="	", header=TRUE, fill=TRUE, comment.char = "")
        universe$Entry.Name<-as.character(universe$Entry.Name)
        

        study<-read.table("/data/_LHG_vs_LCG/(1)LHGvs(-1)LCGincrease.txt", sep="	", header=TRUE, fill=TRUE, comment.char = "")
        study$Entry.Name<-as.character(study$Entry.Name)
        
gf<-GOFrame(GODB, organism="0.01")
gaf<-GOAllFrame(gf)
gsc<-GeneSetCollection(gaf,setType=GOCollection())
params<-GSEAGOHyperGParams(
                name="My params",
                geneSetCollection=gsc,
                geneIds=study$Entry.Name,
                universeGeneIds=universe$Entry.Name,
                pvalueCutoff=1,
                testDirection="over",
                ontology="BP",
                conditional=FALSE)
Over<-hyperGTest(params)
geneIdsByCategory(Over)
summary(Over)
params<-GSEAGOHyperGParams(
                name="My params",
                geneSetCollection=gsc,
                geneIds=study$Entry.Name,
                universeGeneIds=universe$Entry.Name,
                pvalueCutoff=1,
                testDirection="over",
                ontology="CC",
                conditional=FALSE)
Over<-hyperGTest(params)
geneIdsByCategory(Over)
summary(Over)
params<-GSEAGOHyperGParams(
                name="My params",
                geneSetCollection=gsc,
                geneIds=study$Entry.Name,
                universeGeneIds=universe$Entry.Name,
                pvalueCutoff=1,
                testDirection="over",
                ontology="MF",
                conditional=FALSE)
Over<-hyperGTest(params)
geneIdsByCategory(Over)
summary(Over)

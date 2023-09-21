library(GOstats)
library(GSEABase)
GODB<-read.table("/data/_MLTHG_vs_MLTCG/(1)MLTHGvs(-1)MLTCGgostats_association.txt", sep="	", header=TRUE, fill=TRUE, comment.char = "")

        universe<-read.table("/data/_MLTHG_vs_MLTCG/_universe.txt", sep="	", header=TRUE, fill=TRUE, comment.char = "")
        universe$Entry.Name<-as.character(universe$Entry.Name)
        

        study<-read.table("/data/_MLTHG_vs_MLTCG/(1)MLTHGvs(-1)MLTCGdecrease.txt", sep="	", header=TRUE, fill=TRUE, comment.char = "")
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

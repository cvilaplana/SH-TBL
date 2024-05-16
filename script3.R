library(qusage)
library(gdata)

mydata0<-readRDS("FilteredNormalised.rds")

mydata<-t(mydata0)

ph<-read.xls("phenotype.xls")
labs<-ph[,1]


sp<-which(ph$SG !="")





db<-ph[,-1][sp,]

pairs<-db$sample
class<-db$SG

glob<-rep("X",length(pairs))

gcomp<-cbind(glob,db[,c(8:10)])

colnames(gcomp)<-c("G","C","I","E")
ModuleGeneSets = read.gmt("GranulomeModules.gmt")


ModuleNames<-names(ModuleGeneSets)


valuesfc <-matrix(nrow=length(ModuleNames),ncol=dim(gcomp)[2])
valuesfdr<-matrix(1,nrow=length(ModuleNames),ncol=dim(gcomp)[2])

library(ggcorrplot)



for( i in 1:dim(gcomp)[2]) {

SX<-which(gcomp[,i] =="X")

mypairs<-pairs[SX]
myclass<-class[SX]
myset <-mydata[,SX]


res<-qusage(myset,myclass,pairVector=mypairs,"G-H",ModuleGeneSets)

mytable <- qsTable(res)

results<-mytable[match(ModuleNames,mytable[,1]),]

valuesfc[,i]<-results[,2]
valuesfdr[,i]<-results[,4]

}



colnames(valuesfc)<-paste0("log2FC_",colnames(gcomp))
rownames(valuesfc)<-ModuleNames
colnames(valuesfdr)<-paste0("FDR_",colnames(gcomp))
rownames(valuesfdr)<-ModuleNames


Stats<-data.frame(Modules=ModuleNames,valuesfc,valuesfdr)
library(WriteXLS)


WriteXLS("Stats","QusageModuleEnrichment.xls")
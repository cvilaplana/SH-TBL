library(gdata)
library(pheatmap)
library(WGCNA)
library(genefilter)
library(WriteXLS)
library(qusage)

# reading full count matrix 

mydata<-readRDS("TBrnaseq.rds")

mat0<-as.matrix(mydata[,3:dim(mydata)[2]])
rownames(mat0)<-mydata[,2]

mat0<-round(mat0,0)

## from all samples included in the expression matrix we are interested in samples marked as G or H in SG column


ph<-read.xls("phenotype.xls")


sp<-which(ph$SG!="")

mat<-round(mat0[,sp],0)

#db<-ph[,-1][sp,]

#rownames(db)<-colnames(mat)

# same filter as script1

sok<-which(!duplicated(mydata[,2]) & !is.na(mydata[,2]) & (apply(mat,1,max) > 50) & (apply(mat,1,min) > 0) & (apply(mat,1,mean) > 10) )

source("/home/jlozano/JJcode.R")

myset0<-as.matrix(mat[sok,])
rownames(myset0)<-mydata[sok,2]

dfe <- data.frame(groupe = factor(ph$location[sp]))
rownames(dfe)<-colnames(mat)

## Generate a expressed normalised matrix for WGCNA

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = round(myset0, 0), DataFrame(dfe), design = ~groupe)
dds <- DESeq(dds)

# This function calculates a variance stabilizing transformation (VST) from the fitted dispersion-mean relation(s) and then transforms the count data (normalized by division by the size factors or normalization factors), yielding a matrix of values which are now approximately homoskedastic (having constant variance along the range of mean values). The transformation also normalizes with respect to library size.


rld <- vst(dds)
fvals <- assay(rld)


## Generate a filtered by variance normalised matrix for WGCNA taking 10004 genes


myset<-varFilter(fvals, var.func=IQR, var.cutoff=0.3904297, filterByQuantile=TRUE)


exprData<-t(myset)


### We apply WGCNA with a RsquaredCut = 0.8 

powers = c(seq(from = 2, to=30, by=0.5))
sft = pickSoftThreshold(exprData, powerVector = powers, verbose = 5,networkType="signed",RsquaredCut = 0.8 )

power<-sft$powerEstimate # 23 !!

softPower = power;

adjacency = adjacency(exprData, power = softPower, type="signed")

TOM = TOMsimilarity(adjacency,TOMType = "signed")

dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method="average")

minModuleSize=20

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid",deepSplit = 1,minClusterSize = minModuleSize);

dynamicColors = labels2colors(dynamicMods)
nGenes = ncol(exprData);
nSamples = nrow(exprData);

# Recalculate MEs with color labels

MEs0 = moduleEigengenes(exprData,dynamicColors )$eigengenes
MEs = orderMEs(MEs0)


merged<-mergeCloseModules(exprData,dynamicColors,relabel = TRUE)

MEFINAL<-data.frame(merged$newMEs)

#colnames(EMmatrix)<-gsub("MEgrey60","MEorange",colnames(EMmatrix))

#WriteXLS("EMmatrix","EMmatrix.xls")


library(WriteXLS)

# we rename grey60 to orange to discard label problems


merged$colors<-gsub("grey60","orange",merged$colors)
GeneList<-data.frame(colnames(exprData),module=merged$colors)

# WriteXLS("GeneList","GenesAndModules.xls")

### Here we prepare the gmt for Qusage 


um0<-unique(merged$colors)

um<-um0[which(um0 !="grey")]

## 21 modules excluding grey

sink("GranulomeModules.gmt")
for( i in 1:length(um)) {


sxx<-which(merged$colors == um[i] &  um[i] !="grey"  )

cat(um[i],"\tNA\t",paste0(rownames(myset)[sxx],collapse="\t"))

cat("\n")
}

sink()




########################### 

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


for( i in 1:dim(gcomp)[2]) {

SX<-which(gcomp[,i] =="X")

mypairs<-pairs[SX]
myclass<-class[SX]
subset <-myset[,SX]


res<-qusage(subset,myclass,pairVector=mypairs,"G-H",ModuleGeneSets)

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
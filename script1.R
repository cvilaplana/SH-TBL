library(readxl) 
library(pheatmap)
library(DESeq2)
library(WriteXLS)
library(ggpubr)
library(ggbeeswarm)

source("./somefunctions.R")

# reading input files : featurecounts data and some metadata

mydata<-readRDS("./TBrnaseq.rds")
names(mydata) <- gsub("H$", "NL", names(mydata))
mat<-as.matrix(mydata[,3:dim(mydata)[2]])
rownames(mat)<-mydata[,2]
mat <- round(mat, 0)

ph<-as.data.frame(read_excel("./phenotype.xls"))
ph[,2]<-as.factor(ph[,2])

#### (1) Filtering step
# We select unique feature genes and well annotated genes (no NA)  with count in all samples and a counta maximun value of 50 in all samples.

 
sok<-which(!duplicated(mydata[,2]) & !is.na(mydata[,2]) & (apply(mat,1,max) > 50) & (apply(mat,1,min) > 0) & (apply(mat,1,mean) > 10) )

myset<-as.matrix(mat[sok,])
rownames(myset)<-mydata[sok,2]

# DESeq2 preparation

dfa<-data.frame(ph[,3:7])
rownames(dfa)<-ph$SAMPLE

sc<-which(ph$C_NL =="X" | ph$I_NL =="X" | ph$E_NL =="X" )

pho<-data.frame(Group=ph$SG[sc],ph[sc,c(3:7,2)])
rownames(pho)<-ph[sc,1]
phox<-data.frame(Group=ph$SG[sc])
rownames(phox)<-ph[sc,1]


results<-seldeseqU(myset[,sc],phox,"NL","TBL")

SUM40<- results$rnkGeneList[1:40,]
myx<-as.matrix(SUM40[,8:dim(SUM40)[2]])
rownames(myx)<-SUM40[,1]

 ann_colors = list(
         Group = c(NL = "#53BECD", TBL =  "#9D9D9C"),
         Location=  c(NL = "#53BECD", E =  "#60BF49", I = "#FF8000", C ="red"))  

pheatmap(myx,col=rbcol,scale="row",border_color=NA,annotation=phox,fontsize_row=8,fontsize_col=7,annotation_colors = ann_colors[1])

bin<-matrix(phox[,1],nrow=length(phox[,1]))
bin[which(phox[,1]=="NL")]<-0
bin[which(phox[,1]=="TBL")]<-1

TBL40<-rownames(myx)[which(SUM40$log2FoldChange > 1)]
# running ssgsea from 
es<-ssgsea(myset[,sc],list(TBL40))

dfc<-data.frame(group=factor(pho[,2]),ES=as.numeric(es))
comp<-list(c("I","E"),c("NL","C"),c("NL","E"),c("NL","I"),c("C","I"),c("C","E")   )
pal<-c("#53BECD","#60BF49","#FF8000","red")
levs<-c("NL","E","I","C")
dfc<-data.frame(group=factor(pho[,2],levels=levs),ES=as.numeric(es))

P1<-ggboxplot(dfc,x="group",y="ES",palette=pal,fill="group") + xlab("Location") + ylab("Enrichment Score") +  stat_compare_means(comparisons = comp,size=5,method="t.test",label = "p.signif"   )  + geom_beeswarm(cex=2.5,size=2)  + theme(legend.position = "none",axis.text.x = element_text(face="bold", color="black", size=14), axis.text.y =element_text(face="bold", color="black", size=10),axis.title=element_text(size=14,face="bold"),title =element_text(size=16, face='bold'))

###############################################3

### Now we focused in Central Lesions and their respective Non-lesional samples  

sC<-which(ph$C_NL =="X")


phC<-data.frame(Group=ph$location,Ind=ph$sample)[sC,]
rownames(phC)<-ph$SAMPLE[sC]

# running paired DESeq2 by means a wrapper function 

resultsCvsNL <- seldeseqP(myset[,sC],phC,"NL","C")

phox2<-data.frame(Location=ph[sc,3])
rownames(phox2)<-ph[sc,1]

xset<-resultsCvsNL$rnkGeneList[1:40,8:dim(resultsCvsNL$rnkGeneList)[2]]

rownames(xset)<-resultsCvsNL$rnkGeneList[1:40,1]

ann_colors = list(
         Group = c(NL = "#53BECD", TBL =  "#9D9D9C"),
         Location=  c(NL = "#53BECD",C ="red"))  


pheatmap(xset,annotation=phox2,col=rbcol,annotation_colors = ann_colors[2],   scale="row",cluster_col=F,cellwidth=10,fontsize_row=8,border_color=NA)


### Now we focused in Internal Lesions and their respective Non-lesional samples  


sI<-which(ph$I_NL =="X")

phI<-data.frame(Group=ph$location,Ind=ph$sample)[sI,]
rownames(phI)<-ph$SAMPLE[sI]


resultsIvsNL <- seldeseqP(myset[,sI],phI,"NL","I")

xset<-resultsIvsNL$rnkGeneList[1:40,8:dim(resultsIvsNL$rnkGeneList)[2]]
rownames(xset)<-resultsIvsNL$rnkGeneList[1:40,1]


ann_colors = list(
         Group = c(NL = "#53BECD", TBL =  "#9D9D9C"),
         Location=  c(NL = "#53BECD", I = "#FF8000"))  

pheatmap(xset,annotation=phox2,annotation_colors = ann_colors[2],col=rbcol,scale="row",fontsize_row=8,cellwidth=10,cluster_col=F,border_color=NA)



### Now we focused in External Lesions and their respective Non-lesional samples  



sE<-which(ph$E_NL =="X")

phE<-data.frame(Group=ph$location,Ind=ph$sample)[sE,]
rownames(phE)<-ph$SAMPLE[sE]

resultsEvsNL<- seldeseqP(myset[,sE],phE,"NL","E")


xset<-resultsEvsNL$rnkGeneList[1:40,8:dim(resultsEvsNL$rnkGeneList)[2]]
rownames(xset)<-resultsEvsNL$rnkGeneList[1:40,1]

ann_colors = list(
         Group = c(NL = "#53BECD", TBL =  "#9D9D9C"),
         Location=  c(NL = "#53BECD", E =  "#60BF49"))  

pheatmap(xset,annotation=phox2,annotation_colors = ann_colors[2],col=rbcol,scale="row",fontsize_row=7,cellwidth=10,cluster_col=F,border_color=NA)

# note: results$rnkGeneList  include all DESeq2 statistic data from Granulome vs Non-lesional
#
#       resultsEvsNL$rnkGeneList include all DESeq2 statistic data from External vs Non-lesional
#       resultsIvsNL$rnkGeneList include all DESeq2 statistic data from Internall vs Non-lesional
#       resultsCvsNL$rnkGeneList include all DESeq2 statistic data from Central vs Non-lesional
#
rnkGeneList <- as.data.frame(results$rnkGeneList)

write.csv(rnkGeneList, file = "./rnkGeneList.csv")

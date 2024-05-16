library(gdata) # you can change by other R packages as readxl p.ex 
library(pheatmap)
library(DESeq2)
library(WriteXLS)
library(GSVA)
library(ggpubr)
library(ggbeeswarm)

source("somefunctions.R")


# reading input files : featurecounts data and some metadata

mydata<-readRDS("TBrnaseq.rds")
mat<-as.matrix(mydata[,3:dim(mydata)[2]])
rownames(mat)<-mydata[,2]
mat<-round(mat,0)

ph<-read.xls("phenotype.xls")
ph[,2]<-as.factor(ph[,2])



#### (1) Filtering step
# We select unique feature genes and well annotated genes (no NA)  with count in all samples and a counta maximun value of 50 in all samples.

 
sok<-which(!duplicated(mydata[,2]) & !is.na(mydata[,2]) & (apply(mat,1,max) > 50) & (apply(mat,1,min) > 0) & (apply(mat,1,mean) > 10) )

myset<-as.matrix(mat[sok,])
rownames(myset)<-mydata[sok,2]

# DESeq2 preparation

dfa<-data.frame(ph[,3:7])
rownames(dfa)<-ph$SAMPLE

sc<-which(ph$C_H =="X" | ph$I_H =="X" | ph$E_H =="X" )

pho<-data.frame(Group=ph$SG[sc],ph[sc,c(3:7,2)])
rownames(pho)<-ph[sc,1]
phox<-data.frame(Group=ph$SG[sc])
rownames(phox)<-ph[sc,1]


results<-seldeseqU(myset[,sc],phox,"H","G")

SUM40<- results$rnkGeneList[1:40,]
myx<-as.matrix(SUM40[,8:dim(SUM40)[2]])
rownames(myx)<-SUM40[,1]


pdf("Fig2b.pdf")

 ann_colors = list(
         Group = c(H = "cyan", G =  "#F1C40F"),
         Location=  c(H = "cyan", E =  "green", I = "#FF8000", C ="red"))  

pheatmap(myx,col=rbcol,scale="row",border_color=NA,annotation=phox,fontsize_row=8,fontsize_col=7,annotation_colors = ann_colors[1])

dev.off()



bin<-matrix(phox[,1],nrow=length(phox[,1]))
bin[which(phox[,1]=="H")]<-0
bin[which(phox[,1]=="G")]<-1

G30<-rownames(myx)[which(SUM40$log2FoldChange > 1)]
es<-gsva(myset[,sc],list(G30),method="ssgsea")

dfc<-data.frame(group=factor(pho[,2]),ES=as.numeric(es))
comp<-list(c("I","E"),c("H","C"),c("H","E"),c("H","I"),c("C","I"),c("C","E")   )
pal<-c("cyan","green","#FF8000","red")
levs<-c("H","E","I","C")
dfc<-data.frame(group=factor(pho[,2],levels=levs),ES=as.numeric(es))

P1<-ggboxplot(dfc,x="group",y="ES",palette=pal,fill="group") + xlab("Location") + ylab("Enrichment Score") +  stat_compare_means(comparisons = comp,size=5,method="t.test",label = "p.signif"   )  + geom_beeswarm(cex=2.5,size=2)  + theme(legend.position = "none",axis.text.x = element_text(face="bold", color="black", size=14), axis.text.y =element_text(face="bold", color="black", size=10),axis.title=element_text(size=14,face="bold"),title =element_text(size=16, face='bold'))

pdf("Fig2a.pdf")
print(P1)
dev.off()

###############################################3

### Now we focused in Central Lesions and their respective Healthy samples  

sC<-which(ph$C_H =="X")


phC<-data.frame(Group=ph$location,Ind=ph$sample)[sC,]
rownames(phC)<-ph$SAMPLE[sC]

# running paired DESeq2 by means a wrapper function 

resultsCvsH <- seldeseqP(myset[,sC],phC,"H","C")

phox2<-data.frame(Location=ph[sc,3])
rownames(phox2)<-ph[sc,1]

xset<-resultsCvsH$rnkGeneList[1:40,8:dim(resultsCvsH$rnkGeneList)[2]]

rownames(xset)<-resultsCvsH$rnkGeneList[1:40,1]

pdf("Fig3a.pdf")

ann_colors = list(
         Group = c(H = "cyan", G =  "red"),
         Location=  c(H = "cyan",C ="red"))  


pheatmap(xset,annotation=phox2,col=rbcol,annotation_colors = ann_colors[2],   scale="row",cluster_col=F,cellwidth=10,fontsize_row=8,border_color=NA)


### Now we focused in Internal Lesions and their respective Healthy samples  


sI<-which(ph$I_H =="X")

phI<-data.frame(Group=ph$location,Ind=ph$sample)[sI,]
rownames(phI)<-ph$SAMPLE[sI]


resultsIvsH <- seldeseqP(myset[,sI],phI,"H","I")

xset<-resultsIvsH$rnkGeneList[1:40,8:dim(resultsIvsH$rnkGeneList)[2]]
rownames(xset)<-resultsIvsH$rnkGeneList[1:40,1]


ann_colors = list(
         Group = c(H = "cyan", G =  "red"),
         Location=  c(H = "cyan", I = "#FF8000"))  

pheatmap(xset,annotation=phox2,annotation_colors = ann_colors[2],col=rbcol,scale="row",fontsize_row=8,cellwidth=10,cluster_col=F,border_color=NA)



### Now we focused in External Lesions and their respective Healthy samples  



sE<-which(ph$E_H =="X")

phE<-data.frame(Group=ph$location,Ind=ph$sample)[sE,]
rownames(phE)<-ph$SAMPLE[sE]

resultsEvsH<- seldeseqP(myset[,sE],phE,"H","E")


xset<-resultsEvsH$rnkGeneList[1:40,8:dim(resultsEvsH$rnkGeneList)[2]]
rownames(xset)<-resultsEvsH$rnkGeneList[1:40,1]

ann_colors = list(
         Group = c(H = "cyan", G =  "red"),
         Location=  c(H = "cyan", E =  "green"))  



pheatmap(xset,annotation=phox2,annotation_colors = ann_colors[2],col=rbcol,scale="row",fontsize_row=7,cellwidth=10,cluster_col=F,border_color=NA)

dev.off()

# note: results$rnkGeneList  include all DESeq2 statistic data from Granulome vs Healthy
#
#       resultsEvsH$rnkGeneList include all DESeq2 statistic data from External vs Healthy
#       resultsIvsH$rnkGeneList include all DESeq2 statistic data from Internall vs Healthy
#       resultsCvsH$rnkGeneList include all DESeq2 statistic data from Central vs Healthy
# 
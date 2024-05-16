




seldeseqU <- function(cnt,ph,G1,G2) {

sset1<-c(which(ph$Group==G1),c(which(ph$Group==G2)))

mycnt<-cnt[,sset1]
dfe<-data.frame(groupe=factor(ph$Group[sset1],levels=c(G1,G2)))
dds <- DESeqDataSetFromMatrix(countData = round(mycnt,0),DataFrame(dfe),design= ~ groupe)
dds <- DESeq(dds)

myresults <- results(dds)
rld <- rlogTransformation(dds)
fvals <- assay(rld)


dbanot<-data.frame(Gene=rownames(myresults))
GeneList<-data.frame(dbanot,myresults,fvals)

metadata<-data.frame(Group=factor(ph$Group[sset1],levels=c(G1,G2)))



o<-order(myresults$pvalue,decreasing=F)

rnkGeneList<-data.frame(dbanot,myresults,fvals)[o,]

return(list(GeneList=GeneList,rnkGeneList=rnkGeneList,   exprs=fvals,metadata=metadata))

}


seldeseqP<- function(cnt,ph,G1,G2) {

sset1<-c(which(ph$Group==G1),c(which(ph$Group==G2)))

mycnt<-cnt[,sset1]
dfe<-data.frame(groupe=factor(ph$Group[sset1],levels=c(G1,G2)),indv=as.factor(ph$Ind[sset1]))
dds <- DESeqDataSetFromMatrix(countData = round(mycnt,0),DataFrame(dfe),design= ~ indv + groupe)
dds <- DESeq(dds)

myresults <- results(dds)
rld <- rlogTransformation(dds)
fvals <- assay(rld)

GeneList<-data.frame(rownames(cnt),myresults,fvals)
metadata<-data.frame(Group=ph$Group,Indiv=ph$Ind)[sset1,]
rownames(metadata)<-colnames(fvals)

o<-order(myresults$pvalue,decreasing=F)

rnkGeneList<-data.frame(rownames(cnt),myresults,fvals)[o,]

return(list(GeneList=GeneList,rnkGeneList=rnkGeneList,   exprs=fvals,metadata=metadata))

}


rbcol <- colorRampPalette(c("Blue","White","Red"))(256)
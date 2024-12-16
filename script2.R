library(readxl)
library(pheatmap)
library(WGCNA)
library(genefilter)
library(WriteXLS)
library(qusage)

# reading full count matrix 

mydata<-readRDS("./TBrnaseq.rds")

mat0<-as.matrix(mydata[,3:dim(mydata)[2]])
rownames(mat0)<-mydata[,2]

mat0<-round(mat0,0)

## from all samples included in the expression matrix we are interested in samples marked as G or H in SG column


ph<-data.frame(read_excel("./phenotype.xls"))


sp<-which(ph$SG!="")

mat<-round(mat0[,sp],0)

#db<-ph[,-1][sp,]

#rownames(db)<-colnames(mat)

# same filter as script1

sok<-which(!duplicated(mydata[,2]) & !is.na(mydata[,2]) & (apply(mat,1,max) > 50) & (apply(mat,1,min) > 0) & (apply(mat,1,mean) > 10) )


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

sink("./GranulomeModules.gmt")
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

colnames(gcomp)<-c("TBL","C","I","E")
ModuleGeneSets = qusage::read.gmt("./GranulomeModules.gmt")
ModuleNames <- names(ModuleGeneSets)


valuesfc <-matrix(nrow=length(ModuleNames),ncol=dim(gcomp)[2])
valuesfdr<-matrix(1,nrow=length(ModuleNames),ncol=dim(gcomp)[2])


for( i in 1:dim(gcomp)[2]) {
  
  SX<-which(gcomp[,i] =="X")
  
  mypairs<-pairs[SX]
  myclass<-class[SX]
  subset <-myset[,SX]
  
  
  res<-qusage(subset,myclass,pairVector=mypairs,"TBL-NL",ModuleGeneSets)
  
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


WriteXLS("Stats","./QusageModuleEnrichment.xls")

# Locating LM22 genes in modules

library(tidyr)
library(ggplot2)

module_genes <- stack(ModuleGeneSets)
module_genes <- as.data.frame(apply(module_genes, 2, function(x) trimws(x)))
names(module_genes) <- c("Gene", "Module")

lm22 <- read.csv("./LM22.csv")
lm22 <- lm22 %>% 
  pivot_longer(cols = "B_cells_naive":"Neutrophils", names_to = "Population", values_to = "DEG")
lm22 <- lm22[lm22$DEG == 1, ]
lm22 <- lm22[, c("X", "Population")]

lm22_modules <- merge(module_genes, lm22, by.x = "Gene", by.y = "X")
lm22_modules$Population <- sub("(([^_]*_){1}[^_]*)_", "\\1\n", lm22_modules$Population)
lm22_modules$Population <- gsub("_", " ",lm22_modules$Population)
lm22_modules$Population <- gsub(".Tregs.", "(Tregs)",lm22_modules$Population)
lm22_modules$Population <- as.factor(lm22_modules$Population)
lm22_modules <- lm22_modules[lm22_modules$Module != "darkred", ]
lm22_modules$Module <- sapply(lm22_modules$Module, FUN = function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
})
ordered_modules <- c("Salmon", "Brown", "Greenyellow", "Midnightblue", "Yellow", "Turquoise", "Magenta", "Orange", "Lightcyan", "Blue", "Red", "Pink", "Purple", "Green", "Black", "Cyan", "Tan")
ordered_annotations <- c("DNA binding", "EMT", "Neutrophil degranulation", "Cell signaling", "Cell cycle", "Adaptive/Humoral", "Extracellular matrix", "IFN/Cytokine signaling", "Cholesterol biosynthesis", "Metabolism", "Miscellaneous", "Oxidative phosphorylation", "Cilium organization", "Innate/PRR", "Ribosomal/metabolic process", "Myeloid activation", "Organelle biosynthesis")
lm22_modules$Module <- factor(lm22_modules$Module, levels = ordered_modules)
pdf("./SuppFig4.pdf", width = 12, height = 8)
ggplot(lm22_modules, aes(x = Module, fill = Module)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  scale_fill_manual(values = c("salmon", "brown", "greenyellow", "midnightblue",
                               "yellow", "turquoise", "magenta",
                               "orange", "lightcyan", "blue", "red", "pink",
                               "purple", "green", "black", "cyan", "tan"),
                    labels = paste(ordered_modules, ordered_annotations, sep = " - ")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        aspect.ratio = 0.5) +
  labs(y = "Number of genes") +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50), limits = c(0, 50)) +
  facet_wrap(~ Population, ncol = 4)
dev.off()

# Obtaining sub-modules

library(clusterProfiler)
library(org.Hs.eg.db)

turquoise_gene_indices <- which(merged$colors == "turquoise")
turquoise_submodules <- dynamicColors[turquoise_gene_indices]
unique_turquoise_submodules <- unique(turquoise_submodules)
turquoise_submodule_genes <- lapply(unique_turquoise_submodules, function(module_color) {
  gene_indices <- turquoise_gene_indices[turquoise_submodules == module_color]
  gene_names <- colnames(exprData)[gene_indices]
  return(gene_names)
})
names(turquoise_submodule_genes) <- unique_turquoise_submodules

turquoise_submodule_list <- list()
GO_results_list_turquoise <- list()
GO_results_df_list_turquoise <- list()
for (i in 1:length(unique_turquoise_submodules)) {
  turquoise_submodule_list[[i]] <- data.frame(gene = turquoise_submodule_genes[[i]])
  names(turquoise_submodule_list)[i] <- paste0("turquoise_submodule_", i)
  GO_results_list_turquoise[[i]] <- enrichGO(
    gene = turquoise_submodule_list[[i]]$gene,
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    ont = "BP"
  )
  names(GO_results_list_turquoise)[i] <- paste0("GO_results_turquoise_submodule_", i)
  GO_results_df_list_turquoise[[i]] <- as.data.frame(GO_results_list_turquoise[[i]])
  names(GO_results_df_list_turquoise)[i] <- paste0("GO_results_turquoise_submodule_", i, "_df")
}
top10_GO_results_turquoise <- do.call(rbind, lapply(GO_results_df_list_turquoise, function(df) head(df, 10)))
top10_GO_results_turquoise$Submodule <- rep(substr(names(GO_results_df_list_turquoise), 32, 32), each = 10)
row.names(top10_GO_results_turquoise) <- NULL

green_gene_indices <- which(merged$colors == "green")
green_submodules <- dynamicColors[green_gene_indices]
unique_green_submodules <- unique(green_submodules)
green_submodule_genes <- lapply(unique_green_submodules, function(module_color) {
  gene_indices <- green_gene_indices[green_submodules == module_color]
  gene_names <- colnames(exprData)[gene_indices]
  return(gene_names)
})
names(green_submodule_genes) <- unique_green_submodules

green_submodule_list <- list()
GO_results_list_green <- list()
GO_results_df_list_green <- list()
for (i in 1:length(unique_green_submodules)) {
  green_submodule_list[[i]] <- data.frame(gene = green_submodule_genes[[i]])
  names(green_submodule_list)[i] <- paste0("green_submodule_", i)
  GO_results_list_green[[i]] <- enrichGO(
    gene = green_submodule_list[[i]]$gene,
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    ont = "BP"
  )
  names(GO_results_list_green)[i] <- paste0("GO_results_green_submodule_", i)
  GO_results_df_list_green[[i]] <- as.data.frame(GO_results_list_green[[i]])
  names(GO_results_df_list_green)[i] <- paste0("GO_results_green_submodule_", i, "_df")
}
top10_GO_results_green <- do.call(rbind, lapply(GO_results_df_list_green, function(df) head(df, 10)))
top10_GO_results_green$Submodule <- rep(substr(names(GO_results_df_list_green), 32, 32), each = 10)
row.names(top10_GO_results_green) <- NULL

# Sub-module QuSAGE

SubmoduleGeneSets <- module_genes
SubmoduleGeneSets <- SubmoduleGeneSets[SubmoduleGeneSets$Module != "turquoise" & SubmoduleGeneSets$Module != "green", ]
SubmoduleGeneSets <- split(SubmoduleGeneSets$Gene, SubmoduleGeneSets$Module)
for (i in 1:length(names(turquoise_submodule_list))) {
  SubmoduleGeneSets[paste0("turquoise_submodule_", i)] <- turquoise_submodule_list[[i]]
}
for (i in 1:length(names(green_submodule_list))) {
  SubmoduleGeneSets[paste0("green_submodule_", i)] <- green_submodule_list[[i]]
}

SubmoduleNames <- names(SubmoduleGeneSets)

valuesfc <-matrix(nrow=length(SubmoduleNames),ncol=dim(gcomp)[2])
valuesfdr<-matrix(1,nrow=length(SubmoduleNames),ncol=dim(gcomp)[2])

for( i in 1:dim(gcomp)[2]) {
  
  SX<-which(gcomp[,i] =="X")
  
  mypairs<-pairs[SX]
  myclass<-class[SX]
  subset <-myset[,SX]
  
  res<-qusage(subset,myclass,pairVector=mypairs,"TBL-NL",SubmoduleGeneSets)
  
  mytable <- qsTable(res)
  
  results<-mytable[match(SubmoduleNames,mytable[,1]),]
  
  valuesfc[,i]<-results[,2]
  valuesfdr[,i]<-results[,4]
  
}

colnames(valuesfc)<-paste0("log2FC_",colnames(gcomp))
rownames(valuesfc)<-SubmoduleNames
colnames(valuesfdr)<-paste0("FDR_",colnames(gcomp))
rownames(valuesfdr)<-SubmoduleNames

Stats_submodules <-data.frame(Modules=SubmoduleNames,valuesfc,valuesfdr)
WriteXLS("Stats_submodules", "./QusageSubmoduleEnrichment.xls")

# Identifying transcription factors and comparisons between surrogates

library(biomaRt)
library(readxl)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(ggbeeswarm)
library(qusage)

rnkGeneList <- read.csv("./rnkGeneList.csv", row.names = 1)

ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")

GO_ids <- getBM(attributes = c('hgnc_symbol','go_id', 'name_1006'),
                filters = 'hgnc_symbol',
                values = rnkGeneList$Gene,
                mart = ensembl_mart)
tf <- GO_ids[GO_ids$go_id == "GO:0003700", ] #GO term for "DNA binding transcription factor activity"
rnkGeneList_tf <- rnkGeneList[rnkGeneList$Gene %in% tf$hgnc_symbol, ]

ph <- read_excel("./phenotype.xls")
ph <- ph[ph$location != "N" & ph$sample != "3", ] # No NL sample for TB-03
ph <- ph[ph$SAMPLE...1 != "SH.TBL05C", ] # No DESeq data for this sample
ph$SEVERITY <- as.factor(ph$SEVERITY)
ph$Sputum...5 <- as.factor(ph$Sputum...5)

ModuleGeneSets = qusage::read.gmt("./GranulomeModules.gmt")
ModuleGeneSets <- stack(ModuleGeneSets)
ModuleGeneSets <- as.data.frame(apply(ModuleGeneSets, 2, function(x) trimws(x)))
names(ModuleGeneSets) <- c("Gene", "Module")

ModuleGeneSets <- ModuleGeneSets[ModuleGeneSets$Module %in% c("salmon", "brown", "greenyellow", "orange", "lightcyan", "cyan"), ] # Modules associated to surrogates

rnkGeneList_tf <- rnkGeneList_tf[rnkGeneList_tf$Gene %in% ModuleGeneSets$Gene, ]
rnkGeneList_tf <- merge(ModuleGeneSets, rnkGeneList_tf, by = "Gene")
row.names(rnkGeneList_tf) <- rnkGeneList_tf$Gene

tf_ph <- as.data.frame(t(rnkGeneList_tf[, 9:ncol(rnkGeneList_tf)]))
tf_ph$sample <- row.names(tf_ph)
tf_ph <- merge(tf_ph, ph[, c("SAMPLE...1", "SEVERITY", "Sputum...5")], by.x = "sample", by.y = "SAMPLE...1")

tf_names <- names(tf_ph)[2:92]
sev_comparisons <- list()
for (i in tf_names) {
  sev_comparisons[[i]] <- compare_means(as.formula(paste(i, "~ SEVERITY")),
                                        data = tf_ph[, c(2:93)],
                                        method = "wilcox.test",
                                        paired = FALSE,
                                        alternative = "two.sided",
                                        p.adjust.method = "BH")
}
sev_comparisons_df <- do.call(rbind, sev_comparisons)
names(sev_comparisons_df) <- gsub(".y.", "Gene", names(sev_comparisons_df))
sev_comparisons_df <- merge(sev_comparisons_df, rnkGeneList_tf, by = "Gene")
sev_comparisons_df <- sev_comparisons_df[sev_comparisons_df$p.adj <= 0.05, ]
sev_comparisons_df <- sev_comparisons_df[, c("Gene", "Module", "log2FoldChange", "padj","p.adj")]
sev_comparisons_df$comparison <- "Severity"
names(sev_comparisons_df) <- c("Gene", "Module", "Log2FC G vs NL", "p adj. G vs NL", "p adj. Surrogate comparison", "Surrogate comparison")

scc_comparisons <- list()
for (i in tf_names) {
  scc_comparisons[[i]] <- compare_means(as.formula(paste(i, "~ Sputum...5")),
                                        data = tf_ph[, c(2:92, 94)],
                                        method = "wilcox.test",
                                        paired = FALSE,
                                        alternative = "two.sided",
                                        p.adjust.method = "BH")
}
scc_comparisons_df <- do.call(rbind, scc_comparisons)
names(scc_comparisons_df) <- gsub(".y.", "Gene", names(scc_comparisons_df))
scc_comparisons_df <- merge(scc_comparisons_df, rnkGeneList_tf, by = "Gene")
scc_comparisons_df <- scc_comparisons_df[scc_comparisons_df$p.adj <= 0.05, ]
scc_comparisons_df <- scc_comparisons_df[, c("Gene", "Module", "log2FoldChange", "padj","p.adj")]
scc_comparisons_df$comparison <- "SCC"
names(scc_comparisons_df) <- c("Gene", "Module", "Log2FC G vs NL", "p adj. G vs NL", "p adj. Surrogate comparison", "Surrogate comparison")

surr_comparisons <- rbind(sev_comparisons_df, scc_comparisons_df)
surr_comparisons <- surr_comparisons[surr_comparisons$`p adj. G vs NL` <= 0.05, ]
write.csv(surr_comparisons, "./SuppTable2.csv")

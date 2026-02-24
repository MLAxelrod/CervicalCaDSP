#### Overall Cohort ####
### Differential gene expression in high risk v low risk epi and SIME ##

library(ggplot2)
library(rstatix)
library(ggpubr)
library(precrec)
library(cutpointr)
library(ggfortify)
library(DESeq2)
library(fgsea)
library(plyr)
library(tidyverse)
library(ComplexHeatmap)
library(ggrepel)
library(gprofiler2)
library(patchwork)
library(multtest)

theme_set(theme_classic())
theme_update(axis.text=element_text(color="black", size=10), 
             axis.title.y = element_text(size=12), 
             axis.title.x = element_text(size=11))

#zscore function
scaleRow <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
  return(x)
}

path <- "C:/Users/marga/Desktop/CCa_GEOMX/plots/"

cols <- c("#F6222E","#3283FE" , "#FEAF16", "#F8A19F","#2ED9FF", "#90AD1C","#B00068", "#DEA0FD",
          "#325A9B" ,"#C4451C", "#1C8356", "#85660D", "#B10DA1","#1CBE4F", "#C075A6")

##Load data
meta <- read.csv("C:/Users/marga/Desktop/CCa_GEOMX/data/metadata.csv")
rownames(meta) <- paste0("X", meta$SampleID) #R doesn't like row/col names starting with numbers
meta$silva.compart <- paste(meta$Silva, meta$Compartment, sep=".")
table(meta$silva.compart)
meta$ptrisk <- paste(meta$Case, meta$Risk, sep=".")
table(meta$ptrisk)

meta$sil.risk <- paste(meta$Silva, meta$Risk, sep=".")
table(meta$sil.risk, meta$Case)

#Reassign Case ID
meta$CaseID <- revalue(as.character(meta$Case), c("8522" = "WU-01", "11911" = "WU-05", "15381"= "WU-02", "18418"="WU-06", "27635" = "WU-07", "34511"="WU-04", "53047"="WU-03"))


data <- read.csv("C:/Users/marga/Desktop/CCa_GEOMX/data/filtered.norm.counts.linear.formated.csv", row.names = 1)
counts <- read.csv("C:/Users/marga/Desktop/CCa_GEOMX/data/raw.counts.linear.formatted.csv", row.names = 1)

all(rownames(meta)==colnames(data)) ##data to metadata matching
all(rownames(meta)==colnames(counts)) ##data to metadata matching


##Log transform
data <- log(data,2)


##Filter on variable genes
apply(data,1,sd)->row.sd
apply(data,1,mean)->row.mean
row.sd/(row.mean+1e-6)->row.cv 
summary(row.cv)
data.filt <- data[row.cv>0.05,] ##can change threshold

apply(counts,1,sd)->row.sd
apply(counts,1,mean)->row.mean
row.sd/(row.mean+1e-6)->row.cv 
summary(row.cv)
counts.filt <- counts[row.cv>1,]  ##can change threshold

##separate into epi and stroma
meta.epi <- meta[meta$Compartment=="epi",]
counts.epi <- counts.filt[,colnames(counts.filt) %in% rownames(meta.epi)]
data.epi <- data.filt[,colnames(data.filt) %in% rownames(meta.epi)]
all(rownames(meta.epi)==colnames(counts.epi))
all(rownames(meta.epi)==colnames(data.epi))

meta.stroma <- meta[meta$Compartment=="stroma",]
counts.stroma <- counts.filt[,colnames(counts.filt) %in% rownames(meta.stroma)]
data.stroma <- data.filt[,colnames(data.filt) %in% rownames(meta.stroma)]
all(rownames(meta.stroma)==colnames(counts.stroma))
all(rownames(meta.stroma)==colnames(data.stroma))





#########################################################
##### DE Seq by risk in EPI ############################
#######################################################

dds <- DESeqDataSetFromMatrix(countData = counts.epi,
                              colData = meta.epi,
                              design= ~ Risk + CaseID) ### plus caseID to account for non-independent observations

dds$Risk <- relevel(dds$Risk, ref = "Low") # use this to set reference group
dds<-DESeq(dds)
resultsNames(dds)

res<-DESeq2::results(dds, name="Risk_High_vs_Low")
res<-lfcShrink(dds, coef="Risk_High_vs_Low", type="apeglm" )
resOrdered <- res[order(res$pvalue),]
summary(res)
result<-as.data.frame(resOrdered)

result$gene <- rownames(result)
ggplot(result, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  geom_vline(xintercept =0, linetype="dashed")+
  geom_text_repel(data=head(result[result$padj<0.1,],100), aes(label=gene), max.overlaps = 15)+
  labs(x="Log2 Fold Change",
       y= "-Log10(adj. p-value)",
       title="Differentially expressed genes in High-risk v Low-risk epithelium")

#ggsave("HighvLowepi_volc.png", path=path)


#### p value adjusting using two-stage BH procedure
tspadj<- mt.rawp2adjp(result$pvalue, proc="TSBH")
result$padjTSBH <- tspadj$adjp[,"TSBH_0.05"]

#write.csv(result, file = "C:/Users/marga/Desktop/CCa_GEOMX/SR_revision/AllHighvsLowrisk_EPI.csv")

sig<-result[result$padjTSBH<0.1,]
sig<-na.omit(sig)

hivlow <- sig
upvLo <- hivlow[hivlow$log2FoldChange>1,]

##Pathway analysis
gostres <- gost(query=upvLo$gene, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)


gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

#select significant terms of interest
terms <- c("GO:0004175", "GO:0005518", "KEGG:04610", "KEGG:04115", "KEGG:04512", "REAC:R-HSA-1474244", "REAC:R-HSA-1474228", "GO:0008233", "GO:0005201", "GO:0050839")

foo5 <- foo4[foo4$term_id %in% terms,]


ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_point(aes(size=intersection_size, color=source))+
  coord_flip()+
  scale_size(range = c(0,7), limits=c(0,25))+
  scale_color_manual(values= c("#109618","#990099", "#3366cc"))+
  #  ylim(1,15)+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in Higher Risk\nEpithelium (whole cohort)",
       color="Term\nSource",
       size="Intersection\nsize")

#ggsave("pathways_upinHighEpi_wholecohort_NEW.png", path=path)




####################################
##### Risk in Stroma DDS ###########
####################################

dds <- DESeqDataSetFromMatrix(countData = counts.stroma,
                              colData = meta.stroma,
                              design= ~ Risk + CaseID)

dds$Risk <- relevel(dds$Risk, ref = "Low") 
dds<-DESeq(dds)
resultsNames(dds)

## High v Low
res<-DESeq2::results(dds, name="Risk_High_vs_Low")
res<-lfcShrink(dds, coef="Risk_High_vs_Low", type="apeglm" )
resOrdered <- res[order(res$pvalue),]
summary(res)
result<-as.data.frame(resOrdered)

result$gene <- rownames(result)
ggplot(result, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  geom_vline(xintercept =0, linetype="dashed")+
  geom_text_repel(data=head(result[result$padj<0.1,],50), aes(label=gene))+
  labs(x="Log2 Fold Change",
       y= "-Log10(adj. p-value)",
       title="Differentially expressed genes in High-risk v Low-risk stroma")

#ggsave("HighvLowStroma_volc.png", path=path)

#### p value adjusting using two-stage BH procedure
tspadj<- mt.rawp2adjp(result$pvalue, proc="TSBH")
result$padjTSBH <- tspadj$adjp[,"TSBH_0.05"]



#write.csv(result, file = "C:/Users/marga/Desktop/CCa_GEOMX/SR_revision/AllHighvsLowriskStroma.csv")

sig<-result[result$padjTSBH<0.1,]
sig<-na.omit(sig)
hivlow.s <- sig


#### up in high-risk STROMA pathways ####

upvLo <- hivlow.s[hivlow.s$log2FoldChange>0.5,]

gostres <- gost(query=upvLo$gene, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)



foo4 <- gostres$result


terms <- c("GO:0004175", "GO:0005201", "REAC:R-HSA-1474244", "REAC:R-HSA-1592389", "GO:0008233", "GO:0061134", "GO:0005539")


foo5 <- foo4[foo4$term_id %in% terms,]


ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_point(aes(size=intersection_size, color=source))+
  coord_flip()+
  scale_size(range = c(0,5), limits=c(0,20))+
  scale_color_manual(values= c("#109618","#990099", "#3366cc"))+
  #  ylim(1,15)+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in Higher Risk\nSIME (whole cohort)",
       color="Term\nSource",
       size="Intersection\nsize")

#ggsave("pathways_upinHighStroma_wholecohort_new.png", path=path)



############### Overall PCA
pcadata<-cbind(t(data), meta)
pca<-prcomp(t(data), scale=TRUE)
pcadata$Case <- as.character(pcadata$Case)

pcadata$Risk <- revalue(pcadata$Risk, c("Normal" = "Benign"))
pcadata$Compartment <- revalue(pcadata$Compartment, c("epi" = "Epithelium", "stroma" ="SIME"))

autoplot(pca, data= pcadata, colour="Compartment", shape="Risk")
#ggsave("PCAbycompartment.png", path = path)




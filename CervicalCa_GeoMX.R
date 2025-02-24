##CervicalCancer
##GeoMX DSP

###need to make this a renv project ; on reopen some of my packages were gone??

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



data <- read.csv("C:/Users/marga/Desktop/CCa_GEOMX/data/filtered.norm.counts.linear.formated.csv", row.names = 1)
counts <- read.csv("C:/Users/marga/Desktop/CCa_GEOMX/data/raw.counts.linear.formatted.csv", row.names = 1)

all(rownames(meta)==colnames(data)) ##data to metadata matching
all(rownames(meta)==colnames(counts)) ##data to metadata matching

#Hallmarks pathways downloaded from database for GSEA
pathways<- gmtPathways("C:/Users/marga/Desktop/CCa_GEOMX/data/h.all.v2024.1.Hs.symbols.gmt")



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

#### PCA ####

#epi
pcadata<-cbind(t(data.epi), meta.epi)
pca<-prcomp(t(data.epi), scale=TRUE)
pcadata$Case <- as.character(pcadata$Case)
autoplot(pca, data= pcadata, colour= 'Case', shape = 'Risk') +
  # geom_polygon(aes(fill=Risk), alpha = 0.2)+
 # scale_color_manual(values = c( "#1CBE4F","#B10DA1","#FEAF16", "#325A9B"))+
  labs(title="PCA by Patient & Risk in Epithelium")

autoplot(pca, data= pcadata, colour="Silva", shape="Risk")

#ggsave("pca_caserisk_epi.png", path=path)


##stroma
pcadata<-cbind(t(data.stroma), meta.stroma)
pca<-prcomp(t(data.stroma), scale=TRUE)
pcadata$Case <- as.character(pcadata$Case)
autoplot(pca, data= pcadata, colour= 'Case', shape='Risk') +
  # geom_polygon(aes(fill=Risk), alpha = 0.2)+
#  scale_color_manual(values = c( "#1CBE4F","#B10DA1","#FEAF16", "#325A9B"))+
  labs(title="PCA by Patient & Risk in Stroma")

autoplot(pca, data= pcadata, colour="Silva", shape="Risk")
#ggsave("pca_caserisk_stroma.png", path=path)





###DESeq####


##### epithelium v stroma ####
#mainly a proof of concept that there are expected differences in these compartments
dds <- DESeqDataSetFromMatrix(countData = counts.filt,
                              colData = meta,
                              design= ~ Compartment)

dds<-DESeq(dds)
resultsNames(dds)
res<-DESeq2::results(dds, name="Compartment_stroma_vs_epi")
res<-lfcShrink(dds, coef="Compartment_stroma_vs_epi", type="apeglm" )
resOrdered <- res[order(res$pvalue),]
summary(res)
result<-as.data.frame(resOrdered)
result$gene <- rownames(result)


#volcano
ggplot(result, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  geom_vline(xintercept =0, linetype="dashed")+
  geom_text_repel(data=head(result[result$padj<0.1,],50), aes(label=gene))+
  labs(x="Log2 Fold Change",
       y= "-Log10(adj. p-value)",
       title="Differentially expressed genes in Stroma vs. Epithelium")

#ggsave("StromavEpi_volc.png", path=path)

##fgsea
sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways in Stroma vs. Epithelium") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))

#ggsave("gsea_stromavepi.png", path=path)


gseagenes<-fgseaResTidy$leadingEdge[1:3]
foo<-unlist(gseagenes)
foo2<-unique(foo)

up <- sig[sig$log2FoldChange>0.5,]
down <- sig[sig$log2FoldChange<0.5,]

#GO & KEGG analysis
gostres <- gost(query=down$gene, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE)

gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

publish_gosttable(gostres, highlight_terms = gostres$result[c(1:2,10,120),],
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = c("GO:0005737", "GO:0005515"), 
                       width = NA, height = NA, filename = NULL )
pp


##### By Silva in Epithelium ####
dds <- DESeqDataSetFromMatrix(countData = counts.epi,
                              colData = meta.epi,
                              design= ~ Silva)


#dds$Silva <- relevel(dds$Silva, ref = "normal") # use this to set control group, defaults to A

dds<-DESeq(dds)
resultsNames(dds)
## Normal v A
res<-DESeq2::results(dds, name="Silva_normal_vs_A")
res<-lfcShrink(dds, coef="Silva_normal_vs_A", type="apeglm" )
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
       title="Differentially expressed genes in Normal vs Silva A epithelium")

#ggsave("NormalvAepi_volc.png", path=path)


sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways in Normal vs. Silva A") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("gsea_normalvA.png", path=path)


## B v A
res<-DESeq2::results(dds, name="Silva_B_vs_A")
res<-lfcShrink(dds, coef="Silva_B_vs_A", type="apeglm" )
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
       title="Differentially expressed genes in Silva B vs A epithelium")

#ggsave("BvAepi_volc.png", path=path)


sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways in Silva B vs. Silva A") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("gsea_BvA.png", path=path)

## C v A
res<-DESeq2::results(dds, name="Silva_C_vs_A")
res<-lfcShrink(dds, coef="Silva_C_vs_A", type="apeglm" )
resOrdered <- res[order(res$pvalue),]
summary(res)
result<-as.data.frame(resOrdered)
sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways in Silva C vs. Silva A") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("gsea_CvA.png", path=path)

result$gene <- rownames(result)

ggplot(result, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  geom_vline(xintercept =0, linetype="dashed")+
  geom_text_repel(data=head(result[result$padj<0.1,],50), aes(label=gene))+
  labs(x="Log2 Fold Change",
       y= "-Log10(adjusted p-value)",
       title="Differentially expressed genes in Silva C vs A epithelium")
  
#ggsave("CvA_EPI_volc.png", path=path)


### heatmaps
gseagenes<-fgseaResTidy$leadingEdge[1:3]
foo<-unlist(gseagenes)
foo2<-unique(foo)

sig<-data.epi[rownames(data.epi) %in% foo2,]
z<-na.omit(scaleRow(sig))
#z<-na.omit(scaleRow(data.epi))

ann <- data.frame(Silva=meta.epi[["Silva"]])

a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B")))

Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=9),
        show_column_names = FALSE, show_row_dend = FALSE)

z<-na.omit(scaleRow(data.filt))

ann <- data.frame(Silva=meta[["Silva"]],
                  Compartment = meta[["Compartment"]],
                  Patient= as.character(meta[["Case"]]))

a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B"),
                                 Compartment = c("stroma"= "blue", "epi"="red")))

Heatmap(z, 
        top_annotation = a,
        show_row_names = FALSE,
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))




#####DESEQ in stroma by Silva#####

dds <- DESeqDataSetFromMatrix(countData = counts.stroma,
                              colData = meta.stroma,
                              design= ~ Silva)


#dds$Silva <- relevel(dds$Silva, ref = "normal") # use this to set control group, defaults to A

dds<-DESeq(dds)
resultsNames(dds)
## Normal v A
res<-DESeq2::results(dds, name="Silva_normal_vs_A")
res<-lfcShrink(dds, coef="Silva_normal_vs_A", type="apeglm" )
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
       title="Differentially expressed genes in Normal vs. Silva A stroma")

#ggsave("noramlvAstroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways in Normal vs. Silva A Stroma") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("gsea_normalvASTROMA.png", path=path)


## B v A
res<-DESeq2::results(dds, name="Silva_B_vs_A")
res<-lfcShrink(dds, coef="Silva_B_vs_A", type="apeglm" )
resOrdered <- res[order(res$pvalue),]
summary(res)
result<-as.data.frame(resOrdered)
result$gene <- rownames(result)

ggplot(result, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  geom_vline(xintercept =0, linetype="dashed")+
  geom_text_repel(data=head(result[result$padj<0.1,],100), aes(label=gene))+
  labs(x="Log2 Fold Change",
       y= "-Log10(adj. p-value)",
       title="Differentially expressed genes in Silva B vs A Stroma")

#ggsave("BvAstroma_volc.png", path=path)



sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways in Silva B vs. Silva A Stroma") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("gsea_BvASTROMA.png", path=path)

## C v A
res<-DESeq2::results(dds, name="Silva_C_vs_A")
res<-lfcShrink(dds, coef="Silva_C_vs_A", type="apeglm" )
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
       title="Differentially expressed genes in Silva C vs A Stroma")

#ggsave("CvAstrona_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways in Silva C vs. Silva A Stroma") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("gsea_CvASTROMA.png", path=path)

result$gene <- rownames(result)

ggplot(result, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  geom_vline(xintercept =0, linetype="dashed")+
  geom_text_repel(data=head(result[result$padj<0.1,],50), aes(label=gene))+
  labs(x="Log2 Fold Change",
       y= "-Log10(adj. p-value)",
       title="Differentially expressed genes in Silva C vs A Stroma")

#ggsave("CvAstroma_volc.png", path=path)




###############################################################
##### DE Seq by risk in EPI ############
#######################################################

dds <- DESeqDataSetFromMatrix(countData = counts.epi,
                              colData = meta.epi,
                              design= ~ Risk)


dds$Risk <- relevel(dds$Risk, ref = "Normal") # use this to set reference group
dds<-DESeq(dds)
resultsNames(dds)

## High v Normal
res<-DESeq2::results(dds, name="Risk_High_vs_Normal")
res<-lfcShrink(dds, coef="Risk_High_vs_Normal", type="apeglm" )
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
       title="Differentially expressed genes in High-risk v Normal epithelium")

#ggsave("HighvNORMepi_volc.png", path=path)


sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways in High-risk v Normal epithelium") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("gsea_HighvNorm.png", path=path)

hivnorm <- sig

## Low v Normal
res<-DESeq2::results(dds, name="Risk_Low_vs_Normal")
res<-lfcShrink(dds, coef="Risk_Low_vs_Normal", type="apeglm" )
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
       title="Differentially expressed genes in Low-risk v Normal epithelium")

#ggsave("LowvNORMepi_volc.png", path=path)


sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways in Loq-risk v Normal epithelium") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("gsea_LowvNorm.png", path=path)

lovnorm <- sig

#####HIGH V LOW

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
  geom_text_repel(data=head(result[result$padj<0.1,],50), aes(label=gene))+
  labs(x="Log2 Fold Change",
       y= "-Log10(adj. p-value)",
       title="Differentially expressed genes in High-risk v Low-risk epithelium")

#ggsave("HighvLowepi_volc.png", path=path)


sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways in High-risk v Low epithelium") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("gsea_HighvLow.png", path=path)

hivlow <- sig


#### up in high pathways ####

upvLo <- hivlow[hivlow$log2FoldChange>1,]
upvN <- hivnorm[hivnorm$log2FoldChange>1,]
upHi <- unique(c(upvLo$gene, upvN$gene))


gostres <- gost(query=upvLo$gene, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)


gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

#foo5 <- foo4[foo4$intersection_size>15,]
terms <- c("GO:0045236", "GO:0061134", "GO:0008236", "GO:0008233", "GO:0005539", "KEGG:04610", "KEGG:04512",
           "KEGG:04510",  "REAC:R-HSA-1474244", "REAC:R-HSA-1280215", "WP:WP558", "WP:WP129",
           "WP:WP619")


publish_gosttable(gostres, highlight_terms = terms,
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size", "intersection"),
                  filename = NULL)

foo4$per <- foo4$intersection_size/foo4$term_size

p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )+ 
  labs(title = "Up in High-Risk vs Low-Risk Epithelium")


#ggsave("Terms_upinhivlowriskepi.png", path=path)


gostres <- gost(query=upvN$gene, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)


gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result


gostres <- gost(query=upHi, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "KEGG"),
                evcodes = TRUE)

#                sources= c("REAC", "GO:BP", "WP", "KEGG", "GO:MF", "GO:CC"),

##evcodes returns gene list but is slower

gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

terms <- c("KEGG:05322","KEGG:04512", "KEGG:04612", "KEGG:05203", "KEGG:04510", "KEGG:04514",
           "KEGG:04110", "KEGG:05165",
           "REAC:R-HSA-69278","REAC:R-HSA-69002", "REAC:R-HSA-1640170", 
           "REAC:R-HSA-1280215", "REAC:R-HSA-168256")


#terms <- c("GO:0006950", "GO:0050839", "GO:0005515", "GO:0007155", "GO:0002376",
#           "KEGG:05322", "GO:0009611", "GO:0008283", "REAC:R-HSA-69278", 
#           "REAC:R-HSA-69002", "REAC:R-HSA-1640170", "GO:1903047", "GO:0006955",
#           "REAC:R-HSA-1280215", "REAC:R-HSA-168256", "KEGG:04512",
#           "KEGG:04612", "KEGG:05203", "KEGG:04510", "KEGG:04514")


publish_gosttable(gostres, highlight_terms = terms,
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )


pp + labs(title = "Up in High-Risk Epithelium")
#ggsave("GO_upinHi_EPI_keggreac.png", path=path)

foo5 <- foo4[foo4$term_id %in% terms,]

ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_col(aes(fill=source))+
  coord_flip()+
  scale_fill_manual(values= c("#dd4477", "#3366cc"))+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in High-Risk Epithelium",
       fill="Term\nSource")

#ggsave("Hirisk_terms_epi_small.png", path=path)

foo6 <- c(foo5$intersection[1:13])
foo7 <- unlist(strsplit(foo6, split=","))
foo8 <- unique(foo7)


upvN.foo <- upvN[upvN$gene %in% foo8,]
upvL.foo <- upvLo[upvLo$gene %in% foo8,]
upvN.foo2 <- arrange(upvN.foo, desc(log2FoldChange))
genes.keep.epi <- unique(c(upvL.foo$gene, upvN.foo2$gene[1:76]))
#genes.keep.epi <- upvLo$gene



h <- data.epi[rownames(data.epi) %in% genes.keep.epi,]
#h <- data.epi[rownames(data.epi) %in% upHi[1:100],]
z<-na.omit(scaleRow(h))
ann <- data.frame(Risk=meta.epi[["Risk"]])
a <- HeatmapAnnotation(df=ann, 
                       col= list(Risk= c("Low"= "#DEA0FD",
                                          "High"="#F6222E", "Normal"="#325A9B")))




ann <- data.frame(Risk=meta.epi[["Risk"]], Silva=meta.epi[["Silva"]], Patient= as.character(meta.epi[["Case"]]))
a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B"),
                                 Risk= c("Low"= "#DEA0FD",
                                         "High"="#F6222E", "Normal"="#325A9B"),
                                 Patient= c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                            "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                            "8522"= "#DEA0FD")))



Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=6),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))


genes.epi.keggreac <- foo8


h <- data.p11.e[rownames(data.p11.e) %in% genes.keep.epi,]
z<-na.omit(scaleRow(h))
ann <- data.frame(Silva=meta.p11.e[["Silva"]])
a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B")))

Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=7),
        show_column_names = FALSE, show_row_dend = TRUE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))




################################
##### Risk in Stroma DDS ###

dds <- DESeqDataSetFromMatrix(countData = counts.stroma,
                              colData = meta.stroma,
                              design= ~ Risk)


dds$Risk <- relevel(dds$Risk, ref = "Normal") # use this to set reference group
dds<-DESeq(dds)
resultsNames(dds)

## High v Normal
res<-DESeq2::results(dds, name="Risk_High_vs_Normal")
res<-lfcShrink(dds, coef="Risk_High_vs_Normal", type="apeglm" )
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
       title="Differentially expressed genes in High-risk v Normal stroma")

#ggsave("HighvNORMstroma_volc.png", path=path)


sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways in High-risk v Normal stroma") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("gsea_HighvNormStroma.png", path=path)

hivnorm.s <- sig

## Low v Normal
res<-DESeq2::results(dds, name="Risk_Low_vs_Normal")
res<-lfcShrink(dds, coef="Risk_Low_vs_Normal", type="apeglm" )
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
       title="Differentially expressed genes in Low-risk v Normal stroma")

#ggsave("LowvNORMstroma_volc.png", path=path)


sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways in Low-risk v Normal stroma") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("gsea_LowvNormStroma.png", path=path)

lovnorm.s <- sig

#####HIGH V LOW
dds$Risk <- relevel(dds$Risk, ref = "Low") # use this to set reference group
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


sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways in High-risk v Low stroma") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("gsea_HighvLowStroma.png", path=path)

hivlow.s <- sig


#### up in high-risk STROMA pathways ####

upvLo <- hivlow.s[hivlow.s$log2FoldChange>1,]
upvN <- hivnorm.s[hivnorm.s$log2FoldChange>1,]
upHi <- unique(c(upvLo$gene, upvN$gene))




gostres <- gost(query=upvLo$gene, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)



gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result


terms <- c("GO:0061134", "GO:0005201", "REAC:R-HSA-1474244", "GO:0005178")





gostres <- gost(query=upHi, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("KEGG", "REAC"),
                evcodes = TRUE)

#                sources= c("REAC", "GO:BP", "WP", "KEGG", "GO:MF", "GO:CC"),

##evcodes returns gene list but is slower

gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

#terms <- c("GO:0002376", "GO:0006955", "GO:0002682", "GO:0045321", "GO:0006954",
#           "REAC:R-HSA-168256", "GO:0001816", "GO:0042110", "KEGG:04145",
#           "KEGG:04612", "REAC:R-HSA-1280215", "GO:0042605")

terms <- c("KEGG:04145", "KEGG:04612", "KEGG:05322", "KEGG:05165",
           "REAC:R-HSA-168256", "REAC:R-HSA-1280215", "REAC:R-HSA-6798695",
           "REAC:R-HSA-877300", "REAC:R-HSA-168249", "REAC:R-HSA-449147",
           "REAC:R-HSA-1236975", "REAC:R-HSA-913531", "REAC:R-HSA-2132295")



t <- publish_gosttable(gostres, highlight_terms = terms,
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )


pp + labs(title = "Up in High-Risk Stroma")
#ggsave("GO_upinHi_Stroma_keggreac.png", path=path)


foo5 <- foo4[foo4$term_id %in% terms,]

ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_col(aes(fill=source))+
  coord_flip()+
  scale_fill_manual(values= c("#dd4477", "#3366cc"))+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in High-Risk Stroma",
       fill="Term\nSource")

#ggsave("HiRiskTerms_Stroma_condense.png", path=path)

foo6 <- c(foo5$intersection[1:16])
foo7 <- unlist(strsplit(foo6, split=","))
foo8 <- unique(foo7)
genes.stroma.keggreac <- foo8

upvN.foo <- upvN[upvN$gene %in% foo8,]
upvL.foo <- upvLo[upvLo$gene %in% foo8,]
upvN.foo2 <- arrange(upvN.foo, desc(log2FoldChange))
genes.keep.s <- unique(c(upvL.foo$gene, upvN.foo2$gene[1:98]))
genes.keep.s <- upvLo$gene


h <- data.stroma[rownames(data.stroma) %in% genes.keep.s,]
z<-na.omit(scaleRow(h))
ann <- data.frame(Risk=meta.stroma[["Risk"]], Silva=meta.stroma[["Silva"]], Patient= as.character(meta.stroma[["Case"]]))
a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B"),
                                 Risk= c("Low"= "#DEA0FD",
                                         "High"="#F6222E", "Normal"="#325A9B"),
                                 Patient= c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                            "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                            "8522"= "#DEA0FD")))


Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=8),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))




h <- data.p11.s[rownames(data.p11.s) %in% genes.keep.s,]
z<-na.omit(scaleRow(h))
ann <- data.frame(Risk=meta.p11.s[["Risk"]])
a <- HeatmapAnnotation(df=ann, 
                       col= list(Risk= c("Low"= "#DEA0FD",
                                         "High"="#F6222E", "Normal"="#325A9B")))


Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=7),
        show_column_names = FALSE, show_row_dend = TRUE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))











##Stroma genes
Stroma.genes <- unique(c(hivlow.s$gene[1:50], hivnorm.s$gene[1:50], lovnorm.s$gene[1:50]))

h <- data.stroma[rownames(data.stroma) %in% Stroma.genes,]
z<-na.omit(scaleRow(h))
ann <- data.frame(Risk=meta.epi[["Risk"]], Silva=meta.epi[["Silva"]], Patient= as.character(meta.epi[["Case"]]))
a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B"),
                                 Risk= c("Low"= "#DEA0FD",
                                         "High"="#F6222E", "Normal"="#325A9B")))


Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=6),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))


##EPI genes
epi.genes <- unique(c(hivlow$gene[1:50], hivnorm$gene[1:50], lovnorm$gene[1:50]))

h <- data.epi[rownames(data.epi) %in% epi.genes,]
z<-na.omit(scaleRow(h))
ann <- data.frame(Risk=meta.epi[["Risk"]], Silva=meta.epi[["Silva"]], Patient= as.character(meta.epi[["Case"]]))
a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B"),
                                 Risk= c("Low"= "#DEA0FD",
                                         "High"="#F6222E", "Normal"="#325A9B")))


Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=6),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))



#############################################################
########### plots some individual genes #####################

all(rownames(meta)==colnames(data)) 

meta$Case <- as.character(meta$Case)

full <- cbind(meta, as.data.frame(t(scaleRow(data)))) #with Z score

#full <- cbind(meta, as.data.frame(t(data))) # no z score

full$ID2 <- NULL
###really just get rid of this in metadata

full.e <- full[full$Compartment=="epi",]

ggplot(full.e, aes(x= sil.risk, y= SERPINA3))+
  geom_boxplot(outlier.shape = NA, aes(color=Case))+
#  geom_point(aes(color=Case))+
#  geom_line(aes(group=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low Risk", "B High Risk", "C"))+
  labs(x=element_blank(),
       color="Case")


stat <- wilcox_test(full.e, KRT6A~sil.risk)

p1 <- ggplot(full.e, aes(x=sil.risk, y=KRT6A))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="KRT6A Z-score",
       color="Case")+
  theme(legend.position = "none")+
  stat_pvalue_manual(stat, y.position = c(2.3,3.7), hide.ns = TRUE, label="p.adj.signif")


stat <- wilcox_test(full.e, TNC~sil.risk)

p2 <- ggplot(full.e, aes(x=sil.risk, y=TNC))+
  geom_boxplot(outlier.shape = NA)+
    geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="TNC Z-score",
       color="Case")+
  theme(legend.position = "none")+
  stat_pvalue_manual(stat, y.position = c(2,2.5), hide.ns = TRUE, label="p.adj.signif")

stat <- wilcox_test(full.e, SERPINA3 ~ sil.risk)

p3 <- ggplot(full.e, aes(x=sil.risk, y=SERPINA3))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="SERPINA3 Z-score",
       color="Case")+
  theme(legend.position = "none")+
  stat_pvalue_manual(stat, hide.ns = TRUE, label= "p.adj.signif", y.position = c(2.9, 3.7 , 4.1, 3.4, 4.3, 4.6))


stat <- wilcox_test(full.e, LAMC2~sil.risk)

p4 <- ggplot(full.e, aes(x=sil.risk, y=LAMC2))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="LAMC2 Z-score",
       color="Case")+
  theme(legend.position = "none")+
  stat_pvalue_manual(stat, y.position = c(2,2.7,3.7), hide.ns = TRUE, label="p.adj.signif")


stat <- wilcox_test(full.e, FN1~sil.risk)
p5 <- ggplot(full.e, aes(x=sil.risk, y=FN1))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="FN1 Z-score",
       color="Case")+
  theme(legend.position = "none")+
  stat_pvalue_manual(stat, y.position = c(2, 1.6, 1.8 ,2.2), hide.ns = TRUE, label="p.adj.signif")

stat <- wilcox_test(full.e, CCND1 ~sil.risk)

p6 <- ggplot(full.e, aes(x=sil.risk, y=CCND1))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="CCND1 Z-score",
       color="Case")+
  theme(legend.position = "none")+
  stat_pvalue_manual(stat, y.position = c(2, 3.8, 2.5, 3.5, 3), hide.ns = TRUE, label="p.adj.signif")



(p1+p2+p3)/(p4+p5+p6)
#ggsave("EPI6sharedgenes.png", path=path)



p1 <- ggplot(full.e, aes(x=sil.risk, y=KRT6A))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="KRT6A",
       color="Case")+
  theme(legend.position = "none")


p2 <- ggplot(full.e, aes(x=sil.risk, y=KRT5))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="KRT5",
       color="Case")+
  theme(legend.position = "none")


p1|p2

#ggsave("CK56.png", path=path)


risk.genes <- c("FN1", "KRT6A", "LAMC2", "TNC")


risk.data <- full.e[,colnames(full.e) %in% risk.genes]

risk.data$score <- rowSums(risk.data)/length(risk.genes)

full.e$score <- risk.data$score


stat <- wilcox_test(full.e, score~sil.risk)%>% add_significance()

ggplot(full.e, aes(x=sil.risk, y=score))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Normal", "A", "B Low Risk", "B High Risk", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="4 Gene Risk Score")+
  stat_pvalue_manual(stat, y.position = c(3,1,2.8,3.3), label= "p.adj.signif", hide.ns = TRUE)

#ggsave("4genescore.png", path=path)

############### plots genes stroma ##################

full.s <- full[full$Compartment=="stroma",]

#stroma.genes ##id'ed by pathways
#"SFRP2" , "MMP9",   "OLR1"    "FN1" #CTHRC1 ###HK3 excluding case 18418
#"IL1RN"    "CHIT1" "HK3" ## not as good of a separation



stat <- wilcox_test(full.s, CTHRC1~sil.risk)%>% add_significance()

p1 <- ggplot(full.s, aes(x=sil.risk, y=CTHRC1))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="CTHRC1 Z-score")+ theme(legend.position = "none")+
  stat_pvalue_manual(stat, y.position = c(3.5,1,3.2,2.7,3.9,4.2), label= "p.adj.signif", hide.ns = TRUE)


ggplot(full.s, aes(x=sil.risk, y=SFRP2))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="SFRP2 Z-score")+ theme(legend.position = "none")


stat <- wilcox_test(full.s, MMP9~sil.risk)%>% add_significance()

p2 <- ggplot(full.s, aes(x=sil.risk, y=MMP9))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="MMP9 Z-score")+ theme(legend.position = "none")+
  stat_pvalue_manual(stat, y.position = c(2.3,3.1), hide.ns = TRUE, label="p.adj.signif")

p4 <- ggplot(full.s, aes(x=sil.risk, y=OLR1))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="OLR1 Z-score")+ theme(legend.position = "none")


p5 <- ggplot(full.s, aes(x=sil.risk, y=FN1))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="FN1 Z-score")+ theme(legend.position = "none")



ggplot(full.s, aes(x=sil.risk, y=C3))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case")+theme(legend.position = "bottom")

ggsave("legend.png", path=path)


stat <- wilcox_test(full.s, CD68~sil.risk)%>% add_significance()

p7 <- ggplot(full.s, aes(x=sil.risk, y=CD68))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="CD68 Z-score")+ theme(legend.position = "none")+
  stat_pvalue_manual(stat, y.position=c(2.6,3.5,2.3,3,3.3,3.7), hide.ns=TRUE, lable="p.adj.signif")


p1|p2|p7

ggsave("selectStromaGenes.png", path=path)

sub1 <- full.s[full.s$Case!="18418",]

p6 <- ggplot(sub1, aes(x=sil.risk, y=HK3))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case", 
       y="HK3 Z-score")+ theme(legend.position = "none")



(p1+p2+p3)/(p4+p5+p6)


#ggsave("stroma_sixgenes.png", path=path)



####################################
df <- full.e


sepa3 <- df %>% 
  group_by(Case, sil.risk) %>% 
  summarize(mean = mean(SERPINA3),
            min = min(SERPINA3),
            max = max(SERPINA3), 
            q1= quantile(SERPINA3, 0.25),
            q3=quantile(SERPINA3, 0.75))

ggplot(sepa3, aes(x= sil.risk, y=mean, color=Case))+
  geom_boxplot(aes(lower = q1, upper = q3, middle = mean, ymin= min, ymax = max), stat="identity")+
#  geom_line(aes(group=Case))+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Normal", "A", "B Low Risk", "B High Risk", "C"))+
  labs(x=element_blank(), y="SERPINA3 Z-score")
  


ggplot(sepa3, aes(x= sil.risk, y=mean, color=Case))+
  geom_line(aes(group=Case))+
  geom_point()+
  scale_color_manual(values = c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                "8522"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Normal", "A", "B Low Risk", "B High Risk", "C"))+
  labs(x=element_blank(), y="SERPINA3 Z-score")


###do fold change on these 

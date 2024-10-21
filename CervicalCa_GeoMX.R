##CervicalCancer
##GeoMX DSP


### overall: Up in high risk epi = cell cycle, adhesion, immune; Up in high risk stroma is more immune focused
#similar within pt11??
#show pathways and heatmaps?


##To do
## within patient comparisons


##Added high risk and low risk to metadata: LVSI in B cases 34511; 27635

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

#ggsave("pca_caserisk_epi.png", path=path)


##stroma
pcadata<-cbind(t(data.stroma), meta.stroma)
pca<-prcomp(t(data.stroma), scale=TRUE)
pcadata$Case <- as.character(pcadata$Case)
autoplot(pca, data= pcadata, colour= 'Case', shape='Risk') +
  # geom_polygon(aes(fill=Risk), alpha = 0.2)+
#  scale_color_manual(values = c( "#1CBE4F","#B10DA1","#FEAF16", "#325A9B"))+
  labs(title="PCA by Patient & Risk in Stroma")

#ggsave("pca_caserisk_stroma.png", path=path)





###DESeq


##epi v stroma 
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


####By Silva in Epi
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




###DESEQ in stroma###

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



######Within Patient Analyses##############

### Start with 11911 and 34511 which have high low and normal
table(meta$Case) ##first two digits of each case # are unique

##Patient 11 (11911)
meta.p11 <- meta[meta$Case==11911,]
counts.p11 <- counts.filt[,colnames(counts.filt) %in% rownames(meta.p11)]
data.p11 <- data.filt[,colnames(data.filt) %in% rownames(meta.p11)]
all(rownames(meta.p11)==colnames(counts.p11))
all(rownames(meta.p11)==colnames(data.p11))

#PCA
pcadata<-cbind(t(data.p11), meta.p11)
pca<-prcomp(t(data.p11), scale=TRUE)
autoplot(pca, data= pcadata, colour= 'Silva', shape = 'Compartment') +
  # geom_polygon(aes(fill=Risk), alpha = 0.2)+
   scale_color_manual(values = c( "#B10DA1","#FEAF16", "#1CBE4F"))+
  labs(title="Patient 11911")

#ggsave("pca_pt11911.png", path=path)

table(meta.p11$Compartment)
table(meta.p11$silva.compart)

#epi v stroma
meta.p11.e <- meta.p11[meta.p11$Compartment=="epi",]
meta.p11.s <- meta.p11[meta.p11$Compartment=="stroma",]
counts.p11.e <- counts.p11[,colnames(counts.p11) %in% rownames(meta.p11.e)]
counts.p11.s <- counts.p11[,colnames(counts.p11) %in% rownames(meta.p11.s)]
data.p11.e <- data.p11[,colnames(data.p11) %in% rownames(meta.p11.e)]
data.p11.s <- data.p11[,colnames(data.p11) %in% rownames(meta.p11.s)]
all(rownames(meta.p11.e)==colnames(counts.p11.e))
all(rownames(meta.p11.s)==colnames(counts.p11.s))
all(rownames(meta.p11.e)==colnames(data.p11.e))
all(rownames(meta.p11.s)==colnames(data.p11.s))


##DESeq pt 11 epi

dds <- DESeqDataSetFromMatrix(countData = counts.p11.e,
                              colData = meta.p11.e,
                              design= ~ Silva)


#dds$Silva <- relevel(dds$Silva, ref = "normal") # use this to set control group

dds<-DESeq(dds)
resultsNames(dds)

## CvsB
res<-DESeq2::results(dds, name="Silva_C_vs_B")
res<-lfcShrink(dds, coef="Silva_C_vs_B", type="apeglm" )
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
       title="Patient 11911, Silva C vs B in Epithelium")

#ggsave("p11_CvsBepi_volc.png", path=path)

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
       title="Pt 11911, Silva C vs B in epithelium") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("p11_gsea_CvBepi.png", path=path)
p11.CvBepi <- sig

##Comparisons to normal
dds$Silva <- relevel(dds$Silva, ref = "normal") # use this to set control group

dds<-DESeq(dds)
resultsNames(dds)

## Cvsnormal
res<-DESeq2::results(dds, name="Silva_C_vs_normal")
res<-lfcShrink(dds, coef="Silva_C_vs_normal", type="apeglm" )
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
       title="Patient 11911, Silva C vs Normal in Epithelium")

#ggsave("p11_CvsNormepi_volc.png", path=path)

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
       title="Pt 11911, Silva C vs Normal in epithelium") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("p11_gsea_CvNORMepi.png", path=path)

p11.CvNormEpi <- sig

###B vs Normal
res<-DESeq2::results(dds, name="Silva_B_vs_normal")
res<-lfcShrink(dds, coef="Silva_B_vs_normal", type="apeglm" )
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
       title="Patient 11911, Silva B vs Normal in Epithelium")

#ggsave("p11_BvsNormepi_volc.png", path=path)

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
       title="Pt 11911, Silva B vs Normal in epithelium") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("p11_gsea_BvNORMepi.png", path=path)

p11.BvNormEpi <- sig


####GO KEGG on PT11911 comparisions and do stroma comparisons 
###Heatmaps within patient of top 50 genes in each comparison???????

up <- sig[sig$log2FoldChange>0.25,]

#GO & KEGG analysis
gostres <- gost(query=sig$gene, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "GO:BP", "WP", "KEGG"))

gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

terms <- c("REAC:R-HSA-1640170", "REAC:R-HSA-69278", "REAC:R-HSA-69306", "REAC:R-HSA-9018519", 
           "GO:0007049", "GO:1903047", "WP:WP179", "WP:WP466", "KEGG:04110", "KEGG:05203")

publish_gosttable(gostres, highlight_terms = terms,
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )


pp + labs(title = "Pt 11911, Silva B vs. Normal in Epithelium")

#ggsave("p11_GO_BvNorm_epi_condense.png", path=path)


p11.genes <- unique(c(p11.CvBepi$gene[1:40], p11.CvNormEpi$gene[1:40], p11.BvNormEpi$gene[1:40]))

#p11.genes <- unique(c(p11.CvBepi$gene[1:100]))
### heatmaps

p11h <- data.p11.e[rownames(data.p11.e) %in% p11.genes,]
z<-na.omit(scaleRow(p11h))

ann <- data.frame(Silva=meta.p11.e[["Silva"]])

a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B")))

Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=8),
        show_column_names = FALSE, show_row_dend = TRUE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))


h <- data.epi[rownames(data.epi) %in% p11.genes,]
z<-na.omit(scaleRow(h))
ann <- data.frame(Risk=meta.epi[["Risk"]], Silva=meta.epi[["Silva"]], Patient= as.character(meta.epi[["Case"]]))
a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B"),
                                 Risk= c("Low"= "#DEA0FD",
                                         "High"="#F6222E", "Normal"="#325A9B"),
                                 Patient= c("11911"="#F6222E", "15381"="#3283FE" , "18418"="#FEAF16", 
                                            "27635"="#C4451C","34511"="#2ED9FF", "53047"="#1C8356",
                                            "8522"= "#DEA0FD")))


cols <- c("#F6222E","#3283FE" , "#FEAF16", "#F8A19F","#2ED9FF", "#90AD1C","#B00068", "#DEA0FD",
          "#325A9B" ,"#C4451C", "#1C8356", "#85660D", "#B10DA1","#1CBE4F", "#C075A6")


Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=7),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))







###GO on UP IN C EPI
upvB <- p11.CvBepi[p11.CvBepi$log2FoldChange>1,]
upvN <- p11.CvNormEpi[p11.CvNormEpi$log2FoldChange>1,]
upC <- unique(c(upvB$gene, upvN$gene))

gostres <- gost(query=upC, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
#                sources= c("REAC", "GO:BP", "WP", "KEGG", "GO:MF", "GO:CC"),
                sources= c("REAC", "KEGG"),
                evcodes = TRUE)

##evcodes returns gene list

gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

terms <- c("GO:0031982", "GO:0070062", "GO:0002376", "REAC:R-HSA-168256", "GO:0001816", "GO:0006955", "REAC:R-HSA-1280215",
           "WP:WP3888", "REAC:R-HSA-877300", "KEGG:04612", "GO:0048518", "GO:0007155", "GO:0005515",
           "GO:0050839","GO:0045296")

#immune
#terms <- c("GO:0002376", "REAC:R-HSA-168256", "GO:0001816", "GO:0006955", "REAC:R-HSA-1280215", "REAC:R-HSA-877300", "KEGG:04612")
 #adhesion          
#terms <- c("GO:0007155", "GO:0005515","GO:0050839","GO:0045296")


###From overall up in high-risk epi
terms <- c("KEGG:05322","KEGG:04512", "KEGG:04612", "KEGG:05203", "KEGG:04510", "KEGG:04514",
           "KEGG:04110", "KEGG:05165",
           "REAC:R-HSA-69278","REAC:R-HSA-69002", "REAC:R-HSA-1640170", 
           "REAC:R-HSA-1280215", "REAC:R-HSA-168256")


publish_gosttable(gostres, highlight_terms = terms,
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )


pp + labs(title = "Pt 11911, Up in Silva C Epithelium")
#ggsave("pt11_GO_upinSilvaCepi.png", path=path)


foo5 <- foo4[foo4$term_id %in% terms,]
foo6 <- c(foo5$intersection[1:15])
foo7 <- unlist(strsplit(foo6, split=","))
foo8 <- unique(foo7)


p11h <- data.p11.e[rownames(data.p11.e) %in% genes.epi.keggreac,]
z<-na.omit(scaleRow(p11h))
ann <- data.frame(Silva=meta.p11.e[["Silva"]])
a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B")))
Heatmap(z, 
        top_annotation = a,
        show_row_names = FALSE,  row_names_gp = gpar(fontsize=4),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))

### immune and adhesion pathways up in C in pt 11 epithelium


##pt 11 stroma
dds <- DESeqDataSetFromMatrix(countData = counts.p11.s,
                              colData = meta.p11.s,
                              design= ~ Silva)


#dds$Silva <- relevel(dds$Silva, ref = "normal") # use this to set control group

dds<-DESeq(dds)
resultsNames(dds)

## CvsB
res<-DESeq2::results(dds, name="Silva_C_vs_B")
res<-lfcShrink(dds, coef="Silva_C_vs_B", type="apeglm" )
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
       title="Patient 11911, Silva C vs B in Stroma")

#ggsave("p11_CvsBstroma_volc.png", path=path)

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
       title="Pt 11911, Silva C vs B in Stroma") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("p11_gsea_CvBstroma.png", path=path)
p11.CvBstroma <- sig


##Comparisons to normal
dds$Silva <- relevel(dds$Silva, ref = "normal") # use this to set control group

dds<-DESeq(dds)
resultsNames(dds)

## Cvsnormal
res<-DESeq2::results(dds, name="Silva_C_vs_normal")
res<-lfcShrink(dds, coef="Silva_C_vs_normal", type="apeglm" )
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
       title="Patient 11911, Silva C vs Normal in Stroma")

#ggsave("p11_CvsNormStroma_volc.png", path=path)

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
       title="Pt 11911, Silva C vs Normal in Stroma") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("p11_gsea_CvNORMstroma.png", path=path)

p11.CvNormStroma <- sig

###B vs Normal
res<-DESeq2::results(dds, name="Silva_B_vs_normal")
res<-lfcShrink(dds, coef="Silva_B_vs_normal", type="apeglm" )
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
       title="Patient 11911, Silva B vs Normal in Stroma")

#ggsave("p11_BvsNormStroma_volc.png", path=path)

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
       title="Pt 11911, Silva B vs Normal in Stroma") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("p11_gsea_BvNORMstroma.png", path=path)

p11.BvNormStroma <- sig


####GO KEGG on PT11911 comparisions and do stroma comparisons 
###Heatmaps within patient of top 50 genes in each comparison???????

p11.genes.stroma <- unique(c(p11.CvBstroma$gene[1:40], p11.CvNormStroma$gene[1:40], p11.BvNormStroma$gene[1:40]))

up <- p11.CvNormStroma[p11.CvNormStroma$log2FoldChange>0.5,]
down <- p11.CvNormStroma[p11.CvNormStroma$log2FoldChange<0,]

#GO & KEGG analysis
gostres <- gost(query=down$gene, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "GO:BP", "WP", "KEGG", "GO:MF", "GO:CC"))


gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

terms <- c("GO:0002376", "GO:0002682", "GO:0001775", "GO:0007155", "GO:0045321", "GO:0042110",
           "REAC:R-HSA-168256", "GO:0005515", "	REAC:R-HSA-168249", "REAC:R-HSA-1280215", 
           "KEGG:04145", "KEGG:04612")

publish_gosttable(gostres, highlight_terms = terms,
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )


pp + labs(title = "Pt 11911, Silva C vs. Normal in Stroma")

#ggsave("p11_GO_CvNorm_stroma_condense.png", path=path)

p11h <- data.p11.s[rownames(data.p11.s) %in% p11.genes.stroma,]
#p11h <- data.p11.s[rownames(data.p11.s) %in% genes.stroma.keggreac,]
z<-na.omit(scaleRow(p11h))

ann <- data.frame(Silva=meta.p11.s[["Silva"]])

a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B")))

Heatmap(z, 
        top_annotation = a,
        show_row_names = FALSE,  row_names_gp = gpar(fontsize=6),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))



h <- data.stroma[rownames(data.stroma) %in% p11.genes.stroma,]
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


cols <- c("#F6222E","#3283FE" , "#FEAF16", "#F8A19F","#2ED9FF", "#90AD1C","#B00068", "#DEA0FD",
          "#325A9B" ,"#C4451C", "#1C8356", "#85660D", "#B10DA1","#1CBE4F", "#C075A6")


Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=7),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))





####tidy code; do second patient with three and see if similar pathways.. do all? and see ifsimilar pathways
##simplify pathway plots ;come up with more informative heatmap genes... based on pathways? or pathway just the heatmap genes like fewer genes into the pathway thing
# Is there a better like cancer pathways type thing than KEGG and GO?? Maybe better Gene sets for GSEA




### Start with 11911 and 34511 which have high low and normal
table(meta$Case) ##first two digits of each case # are unique

##Patient 34
meta.p34 <- meta[meta$Case==34511,]
counts.p34 <- counts.filt[,colnames(counts.filt) %in% rownames(meta.p34)]
data.p34 <- data.filt[,colnames(data.filt) %in% rownames(meta.p34)]
all(rownames(meta.p34)==colnames(counts.p34))
all(rownames(meta.p34)==colnames(data.p34))

#PCA
pcadata<-cbind(t(data.p34), meta.p34)
pca<-prcomp(t(data.p34), scale=TRUE)
autoplot(pca, data= pcadata, colour= 'Silva', shape = 'Compartment') +
  # geom_polygon(aes(fill=Risk), alpha = 0.2)+
  scale_color_manual(values = c("#1CBE4F","#B10DA1", "#325A9B"))+
  labs(title="Patient 34511")
#ggsave("pca_pt34511.png", path=path)

table(meta.p34$Compartment)
table(meta.p34$silva.compart)

#epi v stroma
meta.p34.e <- meta.p34[meta.p34$Compartment=="epi",]
meta.p34.s <- meta.p34[meta.p34$Compartment=="stroma",]
counts.p34.e <- counts.p34[,colnames(counts.p34) %in% rownames(meta.p34.e)]
counts.p34.s <- counts.p34[,colnames(counts.p34) %in% rownames(meta.p34.s)]
data.p34.e <- data.p34[,colnames(data.p34) %in% rownames(meta.p34.e)]
data.p34.s <- data.p34[,colnames(data.p34) %in% rownames(meta.p34.s)]
all(rownames(meta.p34.e)==colnames(counts.p34.e))
all(rownames(meta.p34.s)==colnames(counts.p34.s))
all(rownames(meta.p34.e)==colnames(data.p34.e))
all(rownames(meta.p34.s)==colnames(data.p34.s))


##DESeq pt 34 epi

dds <- DESeqDataSetFromMatrix(countData = counts.p34.e,
                              colData = meta.p34.e,
                              design= ~ Silva)


#dds$Silva <- relevel(dds$Silva, ref = "normal") # use this to set control group

dds<-DESeq(dds)
resultsNames(dds)

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
       title="Patient 34511, Silva B vs A in Epithelium")

#ggsave("p34_BvAepi_volc.png", path=path)

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
       title="Pt 34511, Silva B vs A in epithelium") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("p34_gsea_BvAepi.png", path=path)
p34.BvAepi <- sig

##Comparisons to normal
dds$Silva <- relevel(dds$Silva, ref = "normal") # use this to set control group

dds<-DESeq(dds)
resultsNames(dds)

## Bvsnormal
res<-DESeq2::results(dds, name="Silva_B_vs_normal")
res<-lfcShrink(dds, coef="Silva_B_vs_normal", type="apeglm" )
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
       title="Patient 34511, Silva B vs Normal in Epithelium")

#ggsave("p34_BvsNormepi_volc.png", path=path)

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
       title="Pt 34511, Silva B vs Normal in epithelium") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("p34_gsea_BvNORMepi.png", path=path)

p34.BvNormEpi <- sig

###A vs Normal
res<-DESeq2::results(dds, name="Silva_A_vs_normal")
res<-lfcShrink(dds, coef="Silva_A_vs_normal", type="apeglm" )
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
       title="Patient 34511, Silva A vs Normal in Epithelium")

#ggsave("p34_AvsNormepi_volc.png", path=path)

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
       title="Pt 11911, Silva A vs Normal in epithelium") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))
#ggsave("p34_gsea_AvNORMepi.png", path=path)

p34.AvNormEpi <- sig


####GO KEGG on PT34 comparisions and do stroma comparisons 


p34.genes <- unique(c(p34.BvAepi$gene[1:40], p34.BvNormEpi$gene[1:40], p34.AvNormEpi$gene[1:40]))

#p34.genes <- unique(c(p11.CvBepi$gene[1:100]))
### heatmaps

p34h <- data.p34.e[rownames(data.p34.e) %in% p34.genes,]
z<-na.omit(scaleRow(p34h))

ann <- data.frame(Silva=meta.p34.e[["Silva"]])

a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B")))

Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=8),
        show_column_names = FALSE, show_row_dend = TRUE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))


###GO on UP IN B
upvA <- p34.BvAepi[p34.BvAepi$log2FoldChange>1,]
upvN <- p34.BvNormEpi[p34.BvNormEpi$log2FoldChange>1,]
upB <- unique(c(upvA$gene, upvN$gene))

gostres <- gost(query=upB, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "GO:BP", "WP", "KEGG", "GO:MF", "GO:CC"),
                evcodes = TRUE)

##evcodes returns gene list but is slower

gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

terms <- c("GO:0070062", "GO:0005515", "REAC:R-HSA-69278", "REAC:R-HSA-1640170", 
           "GO:0000278", "GO:0007049")


publish_gosttable(gostres, highlight_terms = terms,
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )


pp + labs(title = "Pt 34511, Up in Silva B Epithelium")
#ggsave("pt34_GO_upinSilvaBepi.png", path=path)


foo5 <- foo4[foo4$term_id %in% terms,]
foo6 <- c(foo5$intersection[1:15])
foo7 <- unlist(strsplit(foo6, split=","))
foo8 <- unique(foo7)

p34h <- data.p34.e[rownames(data.p34.e) %in% foo8,]
z<-na.omit(scaleRow(p34h))
ann <- data.frame(Silva=meta.p34.e[["Silva"]])
a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B")))
Heatmap(z, 
        top_annotation = a,
        show_row_names = FALSE,  row_names_gp = gpar(fontsize=6),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))

### cell cycle genes up in pt 34 cancer v normal


###Overall do high risk v low risk epi v stroma and then focus on some pt11 highlights

###############################################################
##### DE Seq by risk in EPI ############
######################################

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






















#######################
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



h <- data.stroma[rownames(data.stroma) %in% genes.keep.s,]
z<-na.omit(scaleRow(h))
ann <- data.frame(Risk=meta.stroma[["Risk"]], Silva=meta.stroma[["Silva"]], Patient= as.character(meta.stroma[["Case"]]))
a <- HeatmapAnnotation(df=ann, 
                       col= list(Silva= c("A"= "#1CBE4F", "B" ="#B10DA1",
                                          "C"="#FEAF16", "normal"="#325A9B"),
                                 Risk= c("Low"= "#DEA0FD",
                                         "High"="#F6222E", "Normal"="#325A9B")))


Heatmap(z, 
        top_annotation = a,
        show_row_names = FALSE,  row_names_gp = gpar(fontsize=6),
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




meta$sil.risk <- paste(meta$Silva, meta$Risk, sep=".")
table(meta$sil.risk)
#### add a col for detalied SILVA where you include B with LVSI as separarate or can just rename the col above


###################################
### WITHIN PATIENT ANALYSES ########

library(ggvenn)

### Start with 11911 and 34511 and 53047 which have high low and normal
table(meta$sil.risk, meta$Case) ##first two digits of each case # are unique

meta$cancer <- revalue(meta$Risk, replace= c("High"= "Cancer", "Low"="Cancer"))


########### patient 53 (53047)##########
meta.p53 <- meta[meta$Case==53047,]
counts.p53 <- counts.filt[,colnames(counts.filt) %in% rownames(meta.p53)]
data.p53 <- data.filt[,colnames(data.filt) %in% rownames(meta.p53)]
all(rownames(meta.p53)==colnames(counts.p53))
all(rownames(meta.p53)==colnames(data.p53))
meta.p53.e <- meta.p53[meta.p53$Compartment=="epi",]
meta.p53.s <- meta.p53[meta.p53$Compartment=="stroma",]
counts.p53.e <- counts.p53[,colnames(counts.p53) %in% rownames(meta.p53.e)]
counts.p53.s <- counts.p53[,colnames(counts.p53) %in% rownames(meta.p53.s)]
data.p53.e <- data.p53[,colnames(data.p53) %in% rownames(meta.p53.e)]
data.p53.s <- data.p53[,colnames(data.p53) %in% rownames(meta.p53.s)]
all(rownames(meta.p53.e)==colnames(counts.p53.e))
all(rownames(meta.p53.s)==colnames(counts.p53.s))
all(rownames(meta.p53.e)==colnames(data.p53.e))
all(rownames(meta.p53.s)==colnames(data.p53.s))

## DESeq p53 epi ##
dds <- DESeqDataSetFromMatrix(countData = counts.p53.e,
                              colData = meta.p53.e,
                              design= ~ Silva)

dds<-DESeq(dds)
resultsNames(dds)

## CvsA
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
       title="Patient 53047, Silva C vs A in Epithelium")

#ggsave("p53_CvsAepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p53.CvAepi <- sig

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
       title="Patient 53047, Silva C vs Normal in Epithelium")

#ggsave("p53_CvsNormepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)

p53.CvNormEpi <- sig

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
       title="Patient 53047, Silva A vs Normal in Epithelium")

#ggsave("p53_AvsNormepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)

p53.AvNormEpi <- sig

#####p53 STROMA#####
dds <- DESeqDataSetFromMatrix(countData = counts.p53.s,
                              colData = meta.p53.s,
                              design= ~ Silva)

dds<-DESeq(dds)
resultsNames(dds)

## CvsA
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
       title="Patient 53047, Silva C vs A in Stroma")

#ggsave("p53_CvsAstroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p53.CvAstroma <- sig

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
       title="Patient 53047, Silva C vs Normal in Stroma")

#ggsave("p53_CvsNormStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)

p53.CvNormStroma <- sig

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
       title="Patient 53047, Silva A vs Normal in Stroma")

#ggsave("p53_AvsNormStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)

p53.AvNormStroma <- sig


##cancer vs benign epi
dds <- DESeqDataSetFromMatrix(countData = counts.p53.e,
                              colData = meta.p53.e,
                              design= ~ cancer)

dds$cancer <- relevel(dds$cancer, ref = "Normal")

dds<-DESeq(dds)
resultsNames(dds)

res<-DESeq2::results(dds, name="cancer_Cancer_vs_Normal")
res<-lfcShrink(dds, coef="cancer_Cancer_vs_Normal", type="apeglm" )
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
       title="Patient 53047, Cancer vs Normal in Epithelium")

#ggsave("p53_CanvNormepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p53.CavNormEPI <- sig

##cancer vs benign stroma
dds <- DESeqDataSetFromMatrix(countData = counts.p53.s,
                              colData = meta.p53.s,
                              design= ~ cancer)

dds$cancer <- relevel(dds$cancer, ref = "Normal")

dds<-DESeq(dds)
resultsNames(dds)

res<-DESeq2::results(dds, name="cancer_Cancer_vs_Normal")
res<-lfcShrink(dds, coef="cancer_Cancer_vs_Normal", type="apeglm" )
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
       title="Patient 53047, Cancer vs Normal in Stroma")

#ggsave("p53_CanvNormStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p53.CavNormStroma <- sig

################################



###########Patient 11 (11911)###############
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

p11.BvNormEpi <- sig

#GO & KEGG analysis

p11.upCvBepi <- p11.CvBepi[p11.CvBepi$log2FoldChange>1,]

##Try alternative pathway tools
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

ggo <- groupGO(gene     = p11.upCvBepi$gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE,
               keyType = "SYMBOL")

ego <- enrichGO(gene          = p11.upCvBepi$gene,
#                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE,
                keyType = "SYMBOL")


goplot(ego)+labs(title="pt11 up in C vs B epithelium, GO:MF terms")
#ggsave("p11_CvBepi_GOMF.png", path=path)


###Heatmaps
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


##Apply p11 genes to whole dataset... not sure this analysis makes sense to do
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

Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=7),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))


###GO on UP IN C EPI in p11
upvB <- p11.CvBepi[p11.CvBepi$log2FoldChange>1,]
upvN <- p11.CvNormEpi[p11.CvNormEpi$log2FoldChange>1,]
upC <- unique(c(upvB$gene, upvN$gene))

gostres <- gost(query=upC, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                #                sources= c("REAC", "GO:BP", "WP", "KEGG", "GO:MF", "GO:CC"),
                sources= c("REAC", "KEGG"),
                evcodes = TRUE)

gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

terms <- c("GO:0031982", "GO:0070062", "GO:0002376", "REAC:R-HSA-168256", "GO:0001816", "GO:0006955", "REAC:R-HSA-1280215",
           "WP:WP3888", "REAC:R-HSA-877300", "KEGG:04612", "GO:0048518", "GO:0007155", "GO:0005515",
           "GO:0050839","GO:0045296")

publish_gosttable(gostres, highlight_terms = terms,
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )+ labs(title = "Pt 11911, Up in Silva C Epithelium")
#ggsave("pt11_GO_upinSilvaCepi.png", path=path)


#get genes in interesting terms
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
p11.BvNormStroma <- sig


####GO KEGG on PT11911 stroma
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
                       width = NA, height = NA, filename = NULL )+ labs(title = "Pt 11911, Silva C vs. Normal in Stroma")

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

Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=7),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))


##pt 11 cancer vs normal
dds <- DESeqDataSetFromMatrix(countData = counts.p11.s,
                              colData = meta.p11.s,
                              design= ~ cancer)


dds$cancer <- relevel(dds$cancer, ref = "Normal")

dds<-DESeq(dds)
resultsNames(dds)

res<-DESeq2::results(dds, name="cancer_Cancer_vs_Normal")
res<-lfcShrink(dds, coef="cancer_Cancer_vs_Normal", type="apeglm" )
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
       title="Patient 11911, Cancer vs Normal in Stroma")

#ggsave("p11_CavNormStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p11.CancervNormStroma <- sig

##pt 11 cancer vs normal
dds <- DESeqDataSetFromMatrix(countData = counts.p11.e,
                              colData = meta.p11.e,
                              design= ~ cancer)


dds$cancer <- relevel(dds$cancer, ref = "Normal")

dds<-DESeq(dds)
resultsNames(dds)

res<-DESeq2::results(dds, name="cancer_Cancer_vs_Normal")
res<-lfcShrink(dds, coef="cancer_Cancer_vs_Normal", type="apeglm" )
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
       title="Patient 11911, Cancer vs Normal in Epithelium")

#ggsave("p11_CavNormepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p11.CancervNormepi <- sig


####################################################
############ Patient 34511 ############

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
  geom_text_repel(data=head(result[result$padj<0.1,],200), aes(label=gene))+
  labs(x="Log2 Fold Change",
       y= "-Log10(adj. p-value)",
       title="Patient 34511, Silva B vs A in Epithelium")

#ggsave("p34_BvAepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
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

#####pt 34 stroma###
dds <- DESeqDataSetFromMatrix(countData = counts.p34.s,
                              colData = meta.p34.s,
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
  geom_text_repel(data=head(result[result$padj<0.1,],200), aes(label=gene))+
  labs(x="Log2 Fold Change",
       y= "-Log10(adj. p-value)",
       title="Patient 34511, Silva B vs A in Stroma")

#ggsave("p34_BvAstroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p34.BvAstroma <- sig

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
       title="Patient 34511, Silva B vs Normal in Stroma")

#ggsave("p34_BvsNormStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)

p34.BvNormStroma <- sig

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
       title="Patient 34511, Silva A vs Normal in Stroma")

#ggsave("p34_AvsNormStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)

p34.AvNormStroma <- sig


####pt34 cancer v normal
dds <- DESeqDataSetFromMatrix(countData = counts.p34.e,
                              colData = meta.p34.e,
                              design= ~ cancer)


dds$cancer <- relevel(dds$cancer, ref = "Normal")

dds<-DESeq(dds)
resultsNames(dds)

res<-DESeq2::results(dds, name="cancer_Cancer_vs_Normal")
res<-lfcShrink(dds, coef="cancer_Cancer_vs_Normal", type="apeglm" )
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
       title="Patient 34511, Cancer vs Normal in Epithelium")

#ggsave("p34_CavNormepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p34.CancervNormepi <- sig
##
dds <- DESeqDataSetFromMatrix(countData = counts.p34.s,
                              colData = meta.p34.s,
                              design= ~ cancer)


dds$cancer <- relevel(dds$cancer, ref = "Normal")

dds<-DESeq(dds)
resultsNames(dds)

res<-DESeq2::results(dds, name="cancer_Cancer_vs_Normal")
res<-lfcShrink(dds, coef="cancer_Cancer_vs_Normal", type="apeglm" )
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
       title="Patient 34511, Cancer vs Normal in Stroma")

#ggsave("p34_CavNormStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p34.CancervNormStroma <- sig
##



###look at p11, p53, p34 epi up in high-risk v low risk

up.p11 <- p11.CvBepi[p11.CvBepi$log2FoldChange>1,]
up.p53 <- p53.CvAepi[p53.CvAepi$log2FoldChange>1,]
up.p34 <- p34.BvAepi[p34.BvAepi$log2FoldChange>1,]


library(ggvenn)
a <- list(pt11=up.p11$gene, pt53 = up.p53$gene, pt34 = up.p34$gene)
ggvenn(a) + ggtitle("Shared genes up in high-risk vs low-risk epithelium")
#ggsave("venn_upinhivlowepi.png", path=path)

up.HI <- unique(c(up.p11$gene, up.p34$gene, up.p53$gene))
#710 genes unique
genes.up <- c(up.p11$gene, up.p34$gene, up.p53$gene)
dups <- genes.up[duplicated(genes.up)]
dups2 <- unique(dups)

test <- as.data.frame(table(genes.up))
all <- test[test$Freq>2,]


###STROMA
up.p11.s <- p11.CvBstroma[p11.CvBstroma$log2FoldChange>1,]
up.p53.s <- p53.CvAstroma[p53.CvAstroma$log2FoldChange>1,]
up.p34.s <- p34.BvAstroma[p34.BvAstroma$log2FoldChange>1,]


b <- list(pt11=up.p11.s$gene, pt53 = up.p53.s$gene, pt34 = up.p34.s$gene)

ggvenn(b) + ggtitle("Shared genes up in high-risk vs low-risk stroma")
#ggsave("venn_upinhivlowstroma.png", path=path)


up.HI.s <- unique(c(up.p11.s$gene, up.p34.s$gene, up.p53.s$gene))
genes.up.s <- c(up.p11.s$gene, up.p34.s$gene, up.p53.s$gene)
dups.s <- genes.up.s[duplicated(genes.up.s)]
dups2.s <- unique(dups.s)

test <- as.data.frame(table(genes.up.s))
all <- test[test$Freq>2,] ##No shared genes for stroma all three


########

gostres <- gost(query=dups2.s, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)



gostres <- gost(query=dups2, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)


gostres <- gost(query=up.HI, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)


gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

terms <- c("REAC:R-HSA-8957275", "GO:0005201", "REAC:R-HSA-381426", "GO:0043394", "REAC:R-HSA-1474228", "KEGG:04979", "WP:WP5323")
#above for stroma


terms <- c("REAC:R-HSA-1474244", "GO:0030414", "KEGG:04512", "GO:0050839", "REAC:R-HSA-6785807", "KEGG:05205",
           "GO:0045236", "KEGG:04668", "KEGG:05165", "KEGG:04151")
##above for epi


publish_gosttable(gostres, highlight_terms = terms,
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )+
  labs(title = "Up in High-risk v Low-risk stroma in multiple cases")


#ggsave("terms_upinHivLo_multicases.png", path=path)
#ggsave("terms_upinHivLo_multicases_STROMA.png", path=path)


#write.csv(dups2, file="C:/Users/marga/Desktop/CCa_GEOMX/GenesUpinHivLowEpi_multicase.csv")


h <- data.epi[rownames(data.epi) %in% dups2,]
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

Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=7),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))




##testing other methods on this gene list
library(DOSE)
library(enrichplot)


dups3 <- bitr(dups2, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

edo <- enrichDGN(dups3$ENTREZID)

dotplot(edo, showCategory=20)

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox)
heatplot(edox, showCategory = 5)

edox2 <- pairwise_termsim(edox)
treeplot(edox2)


edo <- pairwise_termsim(edo)
emapplot(edo)

terms <- edo$Description[1:5]
pmcplot(terms, 2010:2020)




#############################################################################
#### More patients #############
#############################################################################

###########Patient 85 (8522)###############
#A, low risk B and normal
meta.p85 <- meta[meta$Case==8522,]
counts.p85 <- counts.filt[,colnames(counts.filt) %in% rownames(meta.p85)]
all(rownames(meta.p85)==colnames(counts.p85))

#epi v stroma
meta.p85.e <- meta.p85[meta.p85$Compartment=="epi",]
meta.p85.s <- meta.p85[meta.p85$Compartment=="stroma",]
counts.p85.e <- counts.p85[,colnames(counts.p85) %in% rownames(meta.p85.e)]
counts.p85.s <- counts.p85[,colnames(counts.p85) %in% rownames(meta.p85.s)]
all(rownames(meta.p85.e)==colnames(counts.p85.e))
all(rownames(meta.p85.s)==colnames(counts.p85.s))


##DESeq pt 85 epi
dds <- DESeqDataSetFromMatrix(countData = counts.p85.e,
                              colData = meta.p85.e,
                              design= ~ Silva)


#dds$Silva <- relevel(dds$Silva, ref = "normal") # use this to set control group

dds<-DESeq(dds)
resultsNames(dds)

## BvsA
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
       title="Patient 8522, Silva B vs A in Epithelium")

#ggsave("p85_BvsAepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p85.BvAepi <- sig

### A v normal
dds$Silva <- relevel(dds$Silva, ref = "normal") # use this to set control group

dds<-DESeq(dds)
resultsNames(dds)
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
       title="Patient 8522, Silva A vs normal in Epithelium")

#ggsave("p85_AvsNormepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p85.Avnormepi <- sig

##### p85 stroma ####

dds <- DESeqDataSetFromMatrix(countData = counts.p85.s,
                              colData = meta.p85.s,
                              design= ~ Silva)

dds<-DESeq(dds)
resultsNames(dds)

## BvsA
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
       title="Patient 8522, Silva B vs A in Stroma")

#ggsave("p85_BvsAstroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p85.BvAstroma <- sig

### A v normal
dds$Silva <- relevel(dds$Silva, ref = "normal") # use this to set control group

dds<-DESeq(dds)
resultsNames(dds)
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
       title="Patient 8522, Silva A vs normal in Stroma")

#ggsave("p85_AvsNormStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p85.AvnormStroma <- sig


####p85 cancer v normal
dds <- DESeqDataSetFromMatrix(countData = counts.p85.e,
                              colData = meta.p85.e,
                              design= ~ cancer)


dds$cancer <- relevel(dds$cancer, ref = "Normal")

dds<-DESeq(dds)
resultsNames(dds)

res<-DESeq2::results(dds, name="cancer_Cancer_vs_Normal")
res<-lfcShrink(dds, coef="cancer_Cancer_vs_Normal", type="apeglm" )
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
       title="Patient 8522, Cancer vs Normal in Epithelium")

#ggsave("p85_CavNormepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p85.CancervNormepi <- sig
##
dds <- DESeqDataSetFromMatrix(countData = counts.p85.s,
                              colData = meta.p85.s,
                              design= ~ cancer)


dds$cancer <- relevel(dds$cancer, ref = "Normal")

dds<-DESeq(dds)
resultsNames(dds)

res<-DESeq2::results(dds, name="cancer_Cancer_vs_Normal")
res<-lfcShrink(dds, coef="cancer_Cancer_vs_Normal", type="apeglm" )
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
       title="Patient 8522, Cancer vs Normal in Stroma")

#ggsave("p85_CavNormStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p85.CancervNormStroma <- sig


##############pt 15381 ##################
#A, low risk B
meta.p15 <- meta[meta$Case==15381,]
counts.p15 <- counts.filt[,colnames(counts.filt) %in% rownames(meta.p15)]
all(rownames(meta.p15)==colnames(counts.p15))

#epi v stroma
meta.p15.e <- meta.p15[meta.p15$Compartment=="epi",]
meta.p15.s <- meta.p15[meta.p15$Compartment=="stroma",]
counts.p15.e <- counts.p15[,colnames(counts.p15) %in% rownames(meta.p15.e)]
counts.p15.s <- counts.p15[,colnames(counts.p15) %in% rownames(meta.p15.s)]
all(rownames(meta.p15.e)==colnames(counts.p15.e))
all(rownames(meta.p15.s)==colnames(counts.p15.s))


##DESeq pt 15 epi
dds <- DESeqDataSetFromMatrix(countData = counts.p15.e,
                              colData = meta.p15.e,
                              design= ~ Silva)

dds<-DESeq(dds)
resultsNames(dds)

## BvsA
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
       title="Patient 15381, Silva B vs A in Epithelium")

#ggsave("p15_BvsAepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p15.BvAepi <- sig

##DESeq pt 15 stroma
dds <- DESeqDataSetFromMatrix(countData = counts.p15.s,
                              colData = meta.p15.s,
                              design= ~ Silva)

dds<-DESeq(dds)
resultsNames(dds)

## BvsA
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
       title="Patient 15381, Silva B vs A in Stroma")

#ggsave("p15_BvsAStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p15.BvAstroma <- sig

#######################pt18418 ################
#C, normal
meta.p18 <- meta[meta$Case==18418,]
counts.p18 <- counts.filt[,colnames(counts.filt) %in% rownames(meta.p18)]
all(rownames(meta.p18)==colnames(counts.p18))

#epi v stroma
meta.p18.e <- meta.p18[meta.p18$Compartment=="epi",]
meta.p18.s <- meta.p18[meta.p18$Compartment=="stroma",]
counts.p18.e <- counts.p18[,colnames(counts.p18) %in% rownames(meta.p18.e)]
counts.p18.s <- counts.p18[,colnames(counts.p18) %in% rownames(meta.p18.s)]
all(rownames(meta.p18.e)==colnames(counts.p18.e))
all(rownames(meta.p18.s)==colnames(counts.p18.s))


##DESeq pt 18 epi
dds <- DESeqDataSetFromMatrix(countData = counts.p18.e,
                              colData = meta.p18.e,
                              design= ~ Silva)


dds$Silva <- relevel(dds$Silva, ref = "normal")
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
       title="Patient 18418, Silva C vs normal in Epithelium")

#ggsave("p18_CvsNormepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p18.Cvnormepi <- sig

##DESeq pt 18 stroma
dds <- DESeqDataSetFromMatrix(countData = counts.p18.s,
                              colData = meta.p18.s,
                              design= ~ Silva)
dds$Silva <- relevel(dds$Silva, ref = "normal")
dds<-DESeq(dds)
resultsNames(dds)

## CvsNorm
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
       title="Patient 18418, Silva C vs normal in Stroma")

#ggsave("p18_CvsNormStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p18.CvNormStroma <- sig


####Cancer vs. Normal epi
dds <- DESeqDataSetFromMatrix(countData = counts.p18.e,
                              colData = meta.p18.e,
                              design= ~ cancer)


dds$cancer <- relevel(dds$cancer, ref = "Normal")

dds<-DESeq(dds)
resultsNames(dds)

res<-DESeq2::results(dds, name="cancer_Cancer_vs_Normal")
res<-lfcShrink(dds, coef="cancer_Cancer_vs_Normal", type="apeglm" )
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
       title="Patient 18418, Cancer vs Normal in Epithelium")

#ggsave("p18_CavNormepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p18.CancervNormepi <- sig


###Cancer v normal stroma
dds <- DESeqDataSetFromMatrix(countData = counts.p18.s,
                              colData = meta.p18.s,
                              design= ~ cancer)


dds$cancer <- relevel(dds$cancer, ref = "Normal")

dds<-DESeq(dds)
resultsNames(dds)

res<-DESeq2::results(dds, name="cancer_Cancer_vs_Normal")
res<-lfcShrink(dds, coef="cancer_Cancer_vs_Normal", type="apeglm" )
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
       title="Patient 18418, Cancer vs Normal in Stroma")

#ggsave("p18_CavNormStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p18.CancervNormStroma <- sig


############pt 27635
#B high vs normal
meta.p27 <- meta[meta$Case==27635,]
counts.p27 <- counts.filt[,colnames(counts.filt) %in% rownames(meta.p27)]
all(rownames(meta.p27)==colnames(counts.p27))

#epi v stroma
meta.p27.e <- meta.p27[meta.p27$Compartment=="epi",]
meta.p27.s <- meta.p27[meta.p27$Compartment=="stroma",]
counts.p27.e <- counts.p27[,colnames(counts.p27) %in% rownames(meta.p27.e)]
counts.p27.s <- counts.p27[,colnames(counts.p27) %in% rownames(meta.p27.s)]
all(rownames(meta.p27.e)==colnames(counts.p27.e))
all(rownames(meta.p27.s)==colnames(counts.p27.s))


##DESeq pt 27 epi
dds <- DESeqDataSetFromMatrix(countData = counts.p27.e,
                              colData = meta.p27.e,
                              design= ~ Silva)


dds$Silva <- relevel(dds$Silva, ref = "normal")
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
       title="Patient 27635, Silva B vs normal in Epithelium")

#ggsave("p27_BvsNormepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p27.Bvnormepi <- sig

##DESeq pt 27 stroma
dds <- DESeqDataSetFromMatrix(countData = counts.p27.s,
                              colData = meta.p27.s,
                              design= ~ Silva)
dds$Silva <- relevel(dds$Silva, ref = "normal")
dds<-DESeq(dds)
resultsNames(dds)

## BvsNorm
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
       title="Patient 27635, Silva B vs normal in Stroma")

#ggsave("p27_BvsNormStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p27.BvNormStroma <- sig

###cancer vs normal epi
dds <- DESeqDataSetFromMatrix(countData = counts.p27.e,
                              colData = meta.p27.e,
                              design= ~ cancer)


dds$cancer <- relevel(dds$cancer, ref = "Normal")

dds<-DESeq(dds)
resultsNames(dds)

res<-DESeq2::results(dds, name="cancer_Cancer_vs_Normal")
res<-lfcShrink(dds, coef="cancer_Cancer_vs_Normal", type="apeglm" )
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
       title="Patient 27635, Cancer vs Normal in Epithelium")

#ggsave("p27_CavNormepi_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p27.CancervNormepi <- sig

#Cancer vs Normal stroma p27
dds <- DESeqDataSetFromMatrix(countData = counts.p27.s,
                              colData = meta.p27.s,
                              design= ~ cancer)


dds$cancer <- relevel(dds$cancer, ref = "Normal")

dds<-DESeq(dds)
resultsNames(dds)

res<-DESeq2::results(dds, name="cancer_Cancer_vs_Normal")
res<-lfcShrink(dds, coef="cancer_Cancer_vs_Normal", type="apeglm" )
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
       title="Patient 27635, Cancer vs Normal in Stroma")

#ggsave("p27_CavNormStroma_volc.png", path=path)

sig<-result[result$padj<0.1,]
sig<-na.omit(sig)
p27.CancervNormStroma <- sig


#################################


########ADD p85 B v normal to do all most invasive vs least invasive 


###Compare all worst pattern of invasion to less invasive
#P85 Blo v A #P11 C vs Blo #P15 Blo v A #P34 Bhi v A #P53 Cv A

#EPI
up.p85 <- p85.BvAepi[p85.BvAepi$log2FoldChange>0.5,]
up.p11 <- p11.CvBepi[p11.CvBepi$log2FoldChange>0.5,]
up.p15 <- p15.BvAepi[p15.BvAepi$log2FoldChange>0.5,]
up.p53 <- p53.CvAepi[p53.CvAepi$log2FoldChange>0.5,]
up.p34 <- p34.BvAepi[p34.BvAepi$log2FoldChange>0.5,]



a <- list(pt11=up.p11$gene, pt53 = up.p53$gene, pt34 = up.p34$gene, pt85=up.p85$gene, pt15=up.p15$gene)
ggvenn(a) + ggtitle("Shared genes up in high-risk vs low-risk epithelium (Log2FC>0.5)")
#ggsave("venn_upinhivlowepi_lfc5.png", path=path)


test <- as.data.frame(table(unlist(a)))

ggplot(test, aes(Freq)) + geom_histogram()+
  geom_text(stat='count', aes(label=..count..), position = position_stack(vjust = 1.2),size=4)+
  labs(x="Number of patients",
       y="Number of genes",
       title="Genes up in worse pattern cancer epithelium, shared across patients")

#ggsave("UpHivLo_epi_histogram.png", path=path)



#### DO go on 3+ vs 2+

up.hi.epi <- test[test$Freq>2,]



gostres <- gost(query=up.hi.epi$Var1, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)


gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

terms <- c("GO:0005201", "GO:0005518", "KEGG:04512", "KEGG:04510", "KEGG:05165", "KEGG:04151", "REAC:R-HSA-1474244", "REAC:R-HSA-1474228", "REAC:R-HSA-3000178")


publish_gosttable(gostres, highlight_terms = terms,
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )+
  labs(title = "Up in High-risk v Low-risk epi in multiple cases")
pp


foo5 <- foo4[foo4$term_id %in% terms,]

c(`GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618", KEGG =
    "#dd4477", REAC = "#3366cc", WP = "#0099c6", TF = "#5574a6", MIRNA = "#22aa99", HPA =
    "#6633cc", CORUM = "#66aa00", HP = "#990099")

ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_col(aes(fill=source))+
  coord_flip()+
 scale_fill_manual(values= c("#109618","#990099", "#3366cc"))+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in Worse Pattern of Invasion in >2 Patients",
       fill="Term\nSource")

#ggsave("Hirisk_terms_epi_small.png", path=path)


foo6 <- foo5[foo5$term_id=="KEGG:04151",]
foo6$intersection

ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_point(aes(size=intersection_size, color=source))+
  coord_flip()+
  scale_size(range = c(0,6), limits=c(0,20))+
  scale_color_manual(values= c("#109618","#990099", "#3366cc"))+
  ylim(1,15)+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in Worse Pattern of Invasion\nEpithelium in 3+ Patients",
       color="Term\nSource",
       size="Intersection\nsize")


#ggsave("terms_upinWorseEPIshared.png", path= path)


#Heatmap

h <- data.epi[rownames(data.epi) %in% up.hi.epi$Var1,]
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

Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=7),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))




########STROMA ################
up.p85 <- p85.BvAstroma[p85.BvAstroma$log2FoldChange>0.5,]
up.p11 <- p11.CvBstroma[p11.CvBstroma$log2FoldChange>0.5,]
up.p15 <- p15.BvAstroma[p15.BvAstroma$log2FoldChange>0.5,]
up.p53 <- p53.CvAstroma[p53.CvAstroma$log2FoldChange>0.5,]
up.p34 <- p34.BvAstroma[p34.BvAstroma$log2FoldChange>0.5,]



a <- list(pt11=up.p11$gene, pt53 = up.p53$gene, pt34 = up.p34$gene, pt85=up.p85$gene, pt15=up.p15$gene)
ggvenn(a) + ggtitle("Shared genes up in high-risk vs low-risk stroma (Log2FC>0.5)")
#ggsave("venn_upinhivlowStroma_lfc5.png", path=path)


test <- as.data.frame(table(unlist(a)))

ggplot(test, aes(Freq)) + geom_histogram()+
  geom_text(stat='count', aes(label=..count..), position = position_stack(vjust = 1.2),size=4)+
  labs(x="Number of patients",
       y="Number of genes",
       title="Genes up in worse pattern cancer stroma,\nshared across patients")

#ggsave("UpHivLo_stroma_histogram.png", path=path)



#### DO go on 3+ vs 2+

up.hi.s <- test[test$Freq>1,]



gostres <- gost(query=up.hi.s$Var1, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)


gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

terms <- c("REAC:R-HSA-6798695", "REAC:R-HSA-168249", "REAC:R-HSA-1474244", "KEGG:04145", "GO:0005125", "WP:WP2806", "GO:0004175", "WP:WP5055")

terms <- c("REAC:R-HSA-6798695", "REAC:R-HSA-168249", "REAC:R-HSA-1474244", "REAC:R-HSA-1474228", "KEGG:04145", "KEGG:04610", "GO:0001968", "GO:0005125", "GO:0004175", "WP:WP2806")


publish_gosttable(gostres, highlight_terms = terms,
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )+
  labs(title = "Up in Worse Pattern Stroma in 2+ cases")
pp


foo5 <- foo4[foo4$term_id %in% terms,]

c(`GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618", KEGG =
    "#dd4477", REAC = "#3366cc", WP = "#0099c6", TF = "#5574a6", MIRNA = "#22aa99", HPA =
    "#6633cc", CORUM = "#66aa00", HP = "#990099")

ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_col(aes(fill=source))+
  coord_flip()+
  scale_fill_manual(values= c("#109618","#990099", "#3366cc", "#dc3912"))+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in Worse Pattern Stroma in 2+ Patients",
       fill="Term\nSource")

#ggsave("WorseStroma_3plusTerms.png", path=path)


ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_point(aes(size=intersection_size, color=source))+
  coord_flip()+
  scale_size(range = c(0,7), limits=c(0,65))+
  scale_color_manual(values= c("#109618","#990099", "#3366cc", "#dc3912"))+
  ylim(1,15)+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in Worse Pattern of Invasion\nStroma in 2+ Patients",
       color="Term\nSource",
       size="Intersection\nsize")


#ggsave("terms_upinWorseStromashared.png", path= path)


###Getting genes driving the upregulated pathways
foo9 <- paste(foo5$intersection, collapse=",")
foo10 <- unlist(strsplit(foo9, split=","))
foo11 <- unique(foo10)

foo12 <- as.character(test$Var1[test$Freq>2])

stroma.genes <- foo11[foo11 %in% foo12]

#Heatmap

h <- data.epi[rownames(data.epi) %in% up.hi.s$Var1,]
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

Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=7),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))










######################## Cancer v Normal epi across pateints ##############################

#EPI
up.p85 <- p85.CancervNormepi[p85.CancervNormepi$log2FoldChange>0.5,]
up.p11 <- p11.CancervNormepi[p11.CancervNormepi$log2FoldChange>0.5,]
up.p53 <- p53.CavNormEPI[p53.CavNormEPI$log2FoldChange>0.5,]
up.p34 <- p34.CancervNormepi[p34.CancervNormepi$log2FoldChange>0.5,]
up.p18 <- p18.CancervNormepi[p18.CancervNormepi$log2FoldChange>0.5,]
up.p27 <- p27.CancervNormepi[p27.CancervNormepi$log2FoldChange>0.5,]

a <- list(pt11=up.p11$gene, pt53 = up.p53$gene, pt34 = up.p34$gene, pt85=up.p85$gene, pt18=up.p18$gene, pt27=up.p27$gene)
ggvenn(a) + ggtitle("Shared genes up in Cancer vs normal")
#ggsave("venn_upinhivlowepi_lfc5.png", path=path)
### NOTE: VENN DOESN'T DRAW FOR MORE THAN 5


test <- as.data.frame(table(unlist(a)))

ggplot(test, aes(Freq)) + geom_histogram()+
  geom_text(stat='count', aes(label=..count..), position = position_stack(vjust = 1.2),size=4)+
  labs(x="Number of patients",
       y="Number of genes",
       title="Genes up in cancer vs normal epithelium,\nshared across patients")

#ggsave("UpCancervNorm_epi_histogram.png", path=path)


# pick shared
up.hi.epi <- test[test$Freq>5,]

gostres <- gost(query=up.hi.epi$Var1, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)


gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

terms <- c("GO:0030527", "GO:0045296", "GO:0050839", "KEGG:04110", "KEGG:03030", "KEGG:05203", "KEGG:04530", "REAC:R-HSA-69278", "REAC:R-HSA-69306",
           "REAC:R-HSA-2559582", "REAC:R-HSA-212300", "WP:WP179", "WP:WP466")


publish_gosttable(gostres, highlight_terms = terms,
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )+
  labs(title = "Up in Cancer vs Normal epi in all cases")
pp


foo5 <- foo4[foo4$term_id %in% terms,]

c(`GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618", KEGG =
    "#dd4477", REAC = "#3366cc", WP = "#0099c6", TF = "#5574a6", MIRNA = "#22aa99", HPA =
    "#6633cc", CORUM = "#66aa00", HP = "#990099")

ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_col(aes(fill=source))+
  coord_flip()+
  scale_fill_manual(values= c("#109618","#990099", "#3366cc", "#dc3912"))+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in Cancer vs Normal in all Patients",
       fill="Term\nSource")

#ggsave("Hirisk_terms_epi_small.png", path=path)


ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_point(aes(size=intersection_size, color=source))+
  coord_flip()+
  scale_size(range = c(0,6), limits=c(0,35))+
  scale_color_manual(values= c("#109618","#990099", "#3366cc", "#dc3912"))+
  ylim(1,22)+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in Cancer vs Normal Epithelium in All Cases",
       color="Term\nSource",
       size="Intersection\nsize")


#ggsave("terms_upinCancervNormEpi_all.png", path= path)


#Heatmap

h <- data.epi[rownames(data.epi) %in% up.hi.epi$Var1,]
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

Heatmap(z, 
        top_annotation = a,
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=7),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))



#####Stroma up in Cancer vs Normal across cases ############
up.p85 <- p85.CancervNormStroma[p85.CancervNormStroma$log2FoldChange>0.5,]
up.p11 <- p11.CancervNormStroma[p11.CancervNormStroma$log2FoldChange>0.5,]
up.p53 <- p53.CavNormStroma[p53.CavNormStroma$log2FoldChange>0.5,]
up.p34 <- p34.CancervNormStroma[p34.CancervNormStroma$log2FoldChange>0.5,]
up.p18 <- p18.CancervNormStroma[p18.CancervNormStroma$log2FoldChange>0.5,]
up.p27 <- p27.CancervNormStroma[p27.CancervNormStroma$log2FoldChange>0.5,]

a <- list(pt11=up.p11$gene, pt53 = up.p53$gene, pt34 = up.p34$gene, pt85=up.p85$gene, pt18=up.p18$gene, pt27=up.p27$gene)
ggvenn(a) + ggtitle("Shared genes up in Cancer vs normal")
#ggsave("venn_upinhivlowepi_lfc5.png", path=path)
### NOTE: VENN DOESN'T DRAW FOR MORE THAN 5


test <- as.data.frame(table(unlist(a)))

ggplot(test, aes(Freq)) + geom_histogram()+
  geom_text(stat='count', aes(label=..count..), position = position_stack(vjust = 1.2),size=4)+
  labs(x="Number of patients",
       y="Number of genes",
       title="Genes up in cancer vs normal stroma,\nshared across patients")

#ggsave("UpCancervNorm_stroma_histogram.png", path=path)


# pick shared
up.hi.s <- test[test$Freq>4,]

gostres <- gost(query=up.hi.s$Var1, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)


gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

terms <- c("GO:0005201", "GO:0042605", "KEGG:04145", "KEGG:04612", "KEGG:04512", "KEGG:04514", "KEGG:05165",
           "REAC:R-HSA-877300", "REAC:R-HSA-1280218", "WP:WP453", "WP:WP306")


publish_gosttable(gostres, highlight_terms = terms,
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )+
  labs(title = "Up in Cancer vs Normal stroma in all cases")
pp


foo5 <- foo4[foo4$term_id %in% terms,]

c(`GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618", KEGG =
    "#dd4477", REAC = "#3366cc", WP = "#0099c6", TF = "#5574a6", MIRNA = "#22aa99", HPA =
    "#6633cc", CORUM = "#66aa00", HP = "#990099")

ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_col(aes(fill=source))+
  coord_flip()+
  scale_fill_manual(values= c("#109618","#990099", "#3366cc", "#dc3912"))+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in Cancer vs Normal in all Patients",
       fill="Term\nSource")

#ggsave("Hirisk_terms_epi_small.png", path=path)


ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_point(aes(size=intersection_size, color=source))+
  coord_flip()+
  scale_size(range = c(0,6), limits=c(0,15))+
  scale_color_manual(values= c("#109618","#990099", "#3366cc", "#dc3912"))+
  ylim(1,11)+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in Cancer vs Normal Stroma in 5+ Cases",
       color="Term\nSource",
       size="Intersection\nsize")


#ggsave("terms_upinCancervNormStroma_shared.png", path= path)


#Heatmap

h <- data.stroma[rownames(data.stroma) %in% up.hi.s$Var1,]
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
        show_row_names = TRUE,  row_names_gp = gpar(fontsize=7),
        show_column_names = FALSE, show_row_dend = FALSE,
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))









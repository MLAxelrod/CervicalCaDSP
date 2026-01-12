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

#p11.upCvBepi <- p11.CvBepi[p11.CvBepi$log2FoldChange>1,]

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




##################################################################################
##############Compare all worst pattern of invasion to less invasive################
#P85 Blo v A #P11 C vs Bhi #P15 Blo v A #P34 Bhi v A #P53 Cv A

#EPI
up.p85 <- p85.BvAepi[p85.BvAepi$log2FoldChange>0.5,]
up.p11 <- p11.CvBepi[p11.CvBepi$log2FoldChange>0.5,]
up.p15 <- p15.BvAepi[p15.BvAepi$log2FoldChange>0.5,]
up.p53 <- p53.CvAepi[p53.CvAepi$log2FoldChange>0.5,]
up.p34 <- p34.BvAepi[p34.BvAepi$log2FoldChange>0.5,]

#write.csv(p85.BvAepi, file="C:/Users/marga/Desktop/CCa_GEOMX/SupplementalData/p85BvAepi.csv")
#write.csv(p11.CvBepi, file="C:/Users/marga/Desktop/CCa_GEOMX/SupplementalData/p11CvBepi.csv")
#write.csv(p15.BvAepi, file="C:/Users/marga/Desktop/CCa_GEOMX/SupplementalData/p15BvAepi.csv")
#write.csv(p53.CvAepi, file="C:/Users/marga/Desktop/CCa_GEOMX/SupplementalData/p53CvAepi.csv")
#write.csv(p34.BvAepi, file="C:/Users/marga/Desktop/CCa_GEOMX/SupplementalData/p34BvAepi.csv")


a <- list(pt11=up.p11$gene, pt53 = up.p53$gene, pt34 = up.p34$gene, pt85=up.p85$gene, pt15=up.p15$gene)
ggvenn(a) + ggtitle("Shared genes up in high-risk vs low-risk epithelium (Log2FC>0.5)")
#ggsave("venn_upinhivlowepi_lfc5.png", path=path)


test <- as.data.frame(table(unlist(a)))

#write.csv(test, file="C:/Users/marga/Desktop/CCa_GEOMX/SupplementalData/GenesUpinHighRiskEpibynumberofcases.csv")

ggplot(test, aes(Freq)) + geom_histogram()+
  geom_text(stat='count', aes(label=..count..), position = position_stack(vjust = 1.2),size=4)+
  labs(x="Number of patients",
       y="Number of genes",
       title="Genes up in worse pattern cancer epithelium, shared across patients")

#ggsave("UpHivLo_epi_histogram.png", path=path)



#### DO go on 3+ vs 2+

up.hi.epi <- test[test$Freq>2,]
#up.hi.epi <- test[test$Freq>1,]



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
       title= "Up in Higher Risk\nEpithelium in 3+ Patients",
       color="Term\nSource",
       size="Intersection\nsize")


#ggsave("terms_upinWorseEPIshared2.png", path= path)



############ up in 2+ pts
up.hi.epi <- test[test$Freq>1,]

gostres <- gost(query=up.hi.epi$Var1, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)


gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

terms <- c("GO:0050839", "GO:0005201", "KEGG:04512", "KEGG:04510", "KEGG:04151", "KEGG:05165", "REAC:R-HSA-1474244", "REAC:R-HSA-1474228", "REAC:R-HSA-216083")


foo5 <- foo4[foo4$term_id %in% terms,]


ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_point(aes(size=intersection_size, color=source))+
  coord_flip()+
  scale_size(range = c(0,6), limits=c(10,75))+
  scale_color_manual(values= c("#109618","#990099", "#3366cc"))+
 ylim(1,35)+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in Higher Risk\nEpithelium in 2+ Patients",
       color="Term\nSource",
       size="Intersection\nsize")


#ggsave("terms_upinWorseEPIshared2plus.png", path= path)





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

# write.csv(p85.BvAstroma, file="C:/Users/marga/Desktop/CCa_GEOMX/SupplementalData/p85BvAstroma.csv")
# write.csv(p11.CvBstroma, file="C:/Users/marga/Desktop/CCa_GEOMX/SupplementalData/p11CvBstroma.csv")
# write.csv(p15.BvAstroma, file="C:/Users/marga/Desktop/CCa_GEOMX/SupplementalData/p15BvAstroma.csv")
# write.csv(p53.CvAstroma, file="C:/Users/marga/Desktop/CCa_GEOMX/SupplementalData/p53CvAstroma.csv")
# write.csv(p34.BvAstroma, file="C:/Users/marga/Desktop/CCa_GEOMX/SupplementalData/p34BvAstroma.csv")



a <- list(pt11=up.p11$gene, pt53 = up.p53$gene, pt34 = up.p34$gene, pt85=up.p85$gene, pt15=up.p15$gene)
ggvenn(a) + ggtitle("Shared genes up in high-risk vs low-risk stroma (Log2FC>0.5)")
#ggsave("venn_upinhivlowStroma_lfc5.png", path=path)


test <- as.data.frame(table(unlist(a)))

#write.csv(test, file="C:/Users/marga/Desktop/CCa_GEOMX/SupplementalData/GenesUpinHighRiskStromabynumberofcases.csv")

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
       title= "Up in Higher Risk Stroma in 2+ Patients",
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
       title= "Up in Higher Risk\nStroma in 2+ Patients",
       color="Term\nSource",
       size="Intersection\nsize")


#ggsave("terms_upinWorseStromashared2.png", path= path)


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









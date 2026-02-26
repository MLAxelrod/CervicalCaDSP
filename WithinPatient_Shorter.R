###################################
### WITHIN PATIENT ANALYSES ########


## Load data and packages from OverallCohort.R
## Code below used to make primary figures

library(ggvenn)

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


#### p value adjusting using two-stage BH procedure
tspadj<- mt.rawp2adjp(result$pvalue, proc="TSBH")
result$padjTSBH <- tspadj$adjp[,"TSBH_0.05"]


sig<-result[result$padjTSBH<0.1,]
sig<-na.omit(sig)
p53.CvAepi <- sig


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

#### p value adjusting using two-stage BH procedure
tspadj<- mt.rawp2adjp(result$pvalue, proc="TSBH")
result$padjTSBH <- tspadj$adjp[,"TSBH_0.05"]


sig<-result[result$padjTSBH<0.1,]
sig<-na.omit(sig)
p53.CvAstroma <- sig



###########Patient 11 (11911)###############
meta.p11 <- meta[meta$Case==11911,]
counts.p11 <- counts.filt[,colnames(counts.filt) %in% rownames(meta.p11)]
data.p11 <- data.filt[,colnames(data.filt) %in% rownames(meta.p11)]
all(rownames(meta.p11)==colnames(counts.p11))
all(rownames(meta.p11)==colnames(data.p11))
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

#### p value adjusting using two-stage BH procedure
tspadj<- mt.rawp2adjp(result$pvalue, proc="TSBH")
result$padjTSBH <- tspadj$adjp[,"TSBH_0.05"]

sig<-result[result$padjTSBH<0.1,]
sig<-na.omit(sig)
p11.CvBepi <- sig

##pt 11 stroma
dds <- DESeqDataSetFromMatrix(countData = counts.p11.s,
                              colData = meta.p11.s,
                              design= ~ Silva)


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

#### p value adjusting using two-stage BH procedure
tspadj<- mt.rawp2adjp(result$pvalue, proc="TSBH")
result$padjTSBH <- tspadj$adjp[,"TSBH_0.05"]

sig<-result[result$padjTSBH<0.1,]
sig<-na.omit(sig)
p11.CvBstroma <- sig


############ Patient 34511 ############
meta.p34 <- meta[meta$Case==34511,]
counts.p34 <- counts.filt[,colnames(counts.filt) %in% rownames(meta.p34)]
data.p34 <- data.filt[,colnames(data.filt) %in% rownames(meta.p34)]
all(rownames(meta.p34)==colnames(counts.p34))
all(rownames(meta.p34)==colnames(data.p34))
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

#### p value adjusting using two-stage BH procedure
tspadj<- mt.rawp2adjp(result$pvalue, proc="TSBH")
result$padjTSBH <- tspadj$adjp[,"TSBH_0.05"]

sig<-result[result$padjTSBH<0.1,]
sig<-na.omit(sig)
p34.BvAepi <- sig

#####pt 34 stroma###
dds <- DESeqDataSetFromMatrix(countData = counts.p34.s,
                              colData = meta.p34.s,
                              design= ~ Silva)

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

#### p value adjusting using two-stage BH procedure
tspadj<- mt.rawp2adjp(result$pvalue, proc="TSBH")
result$padjTSBH <- tspadj$adjp[,"TSBH_0.05"]

sig<-result[result$padjTSBH<0.1,]
sig<-na.omit(sig)
p34.BvAstroma <- sig


###########Patient 85 (8522)###############
#A, low risk B and normal
meta.p85 <- meta[meta$Case==8522,]
counts.p85 <- counts.filt[,colnames(counts.filt) %in% rownames(meta.p85)]
all(rownames(meta.p85)==colnames(counts.p85))
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

#### p value adjusting using two-stage BH procedure
tspadj<- mt.rawp2adjp(result$pvalue, proc="TSBH")
result$padjTSBH <- tspadj$adjp[,"TSBH_0.05"]

sig<-result[result$padjTSBH<0.1,]
sig<-na.omit(sig)
p85.BvAepi <- sig


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

#### p value adjusting using two-stage BH procedure
tspadj<- mt.rawp2adjp(result$pvalue, proc="TSBH")
result$padjTSBH <- tspadj$adjp[,"TSBH_0.05"]

sig<-result[result$padjTSBH<0.1,]
sig<-na.omit(sig)
p85.BvAstroma <- sig


##############pt 15381 ##################
#A, low risk B
meta.p15 <- meta[meta$Case==15381,]
counts.p15 <- counts.filt[,colnames(counts.filt) %in% rownames(meta.p15)]
all(rownames(meta.p15)==colnames(counts.p15))
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

#### p value adjusting using two-stage BH procedure
tspadj<- mt.rawp2adjp(result$pvalue, proc="TSBH")
result$padjTSBH <- tspadj$adjp[,"TSBH_0.05"]

sig<-result[result$padjTSBH<0.1,]
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

#### p value adjusting using two-stage BH procedure
tspadj<- mt.rawp2adjp(result$pvalue, proc="TSBH")
result$padjTSBH <- tspadj$adjp[,"TSBH_0.05"]

sig<-result[result$padjTSBH<0.1,]
sig<-na.omit(sig)
p15.BvAstroma <- sig


##################################################################################
############## Compare all worst pattern of invasion to less invasive ################
#P85 Blo v A #P11 C vs Bhi #P15 Blo v A #P34 Bhi v A #P53 Cv A

#EPITHELIUM 
up.p85 <- p85.BvAepi[p85.BvAepi$log2FoldChange>0.5,]
up.p11 <- p11.CvBepi[p11.CvBepi$log2FoldChange>0.5,]
up.p15 <- p15.BvAepi[p15.BvAepi$log2FoldChange>0.5,]
up.p53 <- p53.CvAepi[p53.CvAepi$log2FoldChange>0.5,]
up.p34 <- p34.BvAepi[p34.BvAepi$log2FoldChange>0.5,]
 
# write.csv(p85.BvAepi, file="C:/Users/marga/Desktop/CCa_GEOMX/SR_revision/p85BvAepi.csv")
# write.csv(p11.CvBepi, file="C:/Users/marga/Desktop/CCa_GEOMX/SR_revision/p11CvBepi.csv")
# write.csv(p15.BvAepi, file="C:/Users/marga/Desktop/CCa_GEOMX/SR_revision/p15BvAepi.csv")
# write.csv(p53.CvAepi, file="C:/Users/marga/Desktop/CCa_GEOMX/SR_revision/p53CvAepi.csv")
# write.csv(p34.BvAepi, file="C:/Users/marga/Desktop/CCa_GEOMX/SR_revision/p34BvAepi.csv")

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


############ Pathway analysis on genes up in 2+ patients higher risk epithelium
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

#ggsave("terms_upinWorseEPIshared2plus_New.png", path= path)



######## STROMA/SIME ################
up.p85 <- p85.BvAstroma[p85.BvAstroma$log2FoldChange>0.5,]
up.p11 <- p11.CvBstroma[p11.CvBstroma$log2FoldChange>0.5,]
up.p15 <- p15.BvAstroma[p15.BvAstroma$log2FoldChange>0.5,]
up.p53 <- p53.CvAstroma[p53.CvAstroma$log2FoldChange>0.5,]
up.p34 <- p34.BvAstroma[p34.BvAstroma$log2FoldChange>0.5,]

# write.csv(p85.BvAstroma, file="C:/Users/marga/Desktop/CCa_GEOMX/SR_revision/p85BvAstroma.csv")
# write.csv(p11.CvBstroma, file="C:/Users/marga/Desktop/CCa_GEOMX/SR_revision/p11CvBstroma.csv")
# write.csv(p15.BvAstroma, file="C:/Users/marga/Desktop/CCa_GEOMX/SR_revision/p15BvAstroma.csv")
# write.csv(p53.CvAstroma, file="C:/Users/marga/Desktop/CCa_GEOMX/SR_revision/p53CvAstroma.csv")
# write.csv(p34.BvAstroma, file="C:/Users/marga/Desktop/CCa_GEOMX/SR_revision/p34BvAstroma.csv")

a <- list(pt11=up.p11$gene, pt53 = up.p53$gene, pt34 = up.p34$gene, pt85=up.p85$gene, pt15=up.p15$gene)
ggvenn(a) + ggtitle("Shared genes up in high-risk vs low-risk stroma (Log2FC>0.5)")
#ggsave("venn_upinhivlowStroma_lfc5.png", path=path)


test <- as.data.frame(table(unlist(a)))

#write.csv(test, file="C:/Users/marga/Desktop/CCa_GEOMX/SR_revision/GenesUpinHighRiskStromabynumberofcases.csv")

ggplot(test, aes(Freq)) + geom_histogram()+
  geom_text(stat='count', aes(label=..count..), position = position_stack(vjust = 1.2),size=4)+
  labs(x="Number of patients",
       y="Number of genes",
       title="Genes up in worse pattern cancer stroma,\nshared across patients")

#ggsave("UpHivLo_stroma_histogram.png", path=path)


#### Pathway analysis on genes up in 2+ cases in higher risk SIME

up.hi.s <- test[test$Freq>1,]

gostres <- gost(query=up.hi.s$Var1, organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                sources= c("REAC", "WP", "KEGG", "GO:MF"),
                evcodes = TRUE)


#gostplot(gostres, capped = FALSE, interactive = TRUE)

foo4 <- gostres$result

terms <- c("REAC:R-HSA-6798695", "REAC:R-HSA-168249", "REAC:R-HSA-1474244", "KEGG:04145", "GO:0005125", "WP:WP2806", "GO:0004175", "WP:WP5055")
terms <- c("REAC:R-HSA-6798695", "REAC:R-HSA-168249", "REAC:R-HSA-1474244", "REAC:R-HSA-1474228", "KEGG:04145", "KEGG:04610", "GO:0001968", "GO:0005125", "GO:0004175", "WP:WP2806")


foo5 <- foo4[foo4$term_id %in% terms,]


ggplot(foo5, aes(reorder(term_name, p_value), -log10(p_value)))+
  geom_point(aes(size=intersection_size, color=source))+
  coord_flip()+
  scale_size(range = c(0,7), limits=c(0,65))+
  scale_color_manual(values= c("#109618","#990099", "#3366cc", "#dc3912"))+
#  ylim(1,15)+
  labs(x="Term Name",
       y="-Log10(p value)",
       title= "Up in Higher Risk\nSIME in 2+ Patients",
       color="Term\nSource",
       size="Intersection\nsize")


#ggsave("terms_upinWorseStromashared2_NEW.png", path= path)



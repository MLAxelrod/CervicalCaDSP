#### Overall Cohort ####
### Differential gene expression in high risk v low risk epi and SIME ##
## Code below used for primary figures #######

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


#### Load data and metadata
meta <- read.csv("C:/Users/marga/Desktop/CCa_GEOMX/data/metadata.csv")
rownames(meta) <- paste0("X", meta$SampleID) #R doesn't like row/col names starting with numbers
meta$silva.compart <- paste(meta$Silva, meta$Compartment, sep=".")
meta$ptrisk <- paste(meta$Case, meta$Risk, sep=".")
meta$sil.risk <- paste(meta$Silva, meta$Risk, sep=".")
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



############### Overall PCA ##################
pcadata<-cbind(t(data), meta)
pca<-prcomp(t(data), scale=TRUE)
pcadata$Case <- as.character(pcadata$Case)
pcadata$Risk <- revalue(pcadata$Risk, c("Normal" = "Benign"))
pcadata$Compartment <- revalue(pcadata$Compartment, c("epi" = "Epithelium", "stroma" ="SIME"))

autoplot(pca, data= pcadata, colour="Compartment", shape="Risk")
#ggsave("PCAbycompartment.png", path = path)


########### Plot individual genes #####################

all(rownames(meta)==colnames(data)) 

meta$Case <- as.character(meta$Case)
meta$ID2 <- NULL

full <- cbind(meta, as.data.frame(t(scaleRow(data)))) #with Z score

#full <- cbind(meta, as.data.frame(t(data))) # no z score


full.e <- full[full$Compartment=="epi",]


stat <- wilcox_test(full.e, KRT6A~sil.risk)

p1 <- ggplot(full.e, aes(x=sil.risk, y=KRT6A))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=CaseID))+
  scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                "WU-07"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="KRT6A Z-score",
       color="Case",
       title="KRT6A")+
  theme(legend.position = "none", plot.title= element_text(face="bold", hjust=0.5))+
  stat_pvalue_manual(stat, y.position = c(2.3,1.8,2.5,3,2.7,3.7), hide.ns = TRUE, label="p.adj.signif")



stat <- wilcox_test(full.e, TNC~sil.risk)

p2 <- ggplot(full.e, aes(x=sil.risk, y=TNC))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=CaseID))+
  scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                "WU-07"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="TNC Z-score",
       color="Case",
       title="TNC")+
  theme(legend.position = "none", plot.title= element_text(face="bold", hjust=0.5))+
  stat_pvalue_manual(stat, y.position = c(1.8,2,2.3), hide.ns = TRUE, label="p.adj.signif")

stat <- wilcox_test(full.e, SERPINA3 ~ sil.risk)

p3 <- ggplot(full.e, aes(x=sil.risk, y=SERPINA3))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=CaseID))+
  scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                "WU-07"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="SERPINA3 Z-score",
       color="Case",
       title="SERPINA3")+
  theme(legend.position = "none", plot.title= element_text(face="bold", hjust=0.5))+
  stat_pvalue_manual(stat, hide.ns = TRUE, label= "p.adj.signif", y.position = c(2.9, 4.3 , 4.6))


stat <- wilcox_test(full.e, LAMC2~sil.risk)

p4 <- ggplot(full.e, aes(x=sil.risk, y=LAMC2))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=CaseID))+
  scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                "WU-07"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="LAMC2 Z-score",
       color="Case",
       title="LAMC2")+
  theme(legend.position = "none", plot.title= element_text(face="bold", hjust=0.5))+
  stat_pvalue_manual(stat, y.position = c(2,2.7,3.7), hide.ns = TRUE, label="p.adj.signif")


stat <- wilcox_test(full.e, FN1~sil.risk)
p5 <- ggplot(full.e, aes(x=sil.risk, y=FN1))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=CaseID))+
  scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                "WU-07"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="FN1 Z-score",
       color="Case",
       title="FN1")+
  theme(legend.position = "none", plot.title= element_text(face="bold", hjust=0.5))+
  stat_pvalue_manual(stat, y.position = c(2, 1.6, 1.8), hide.ns = TRUE, label="p.adj.signif")

stat <- wilcox_test(full.e, CCND1 ~sil.risk)

p6 <- ggplot(full.e, aes(x=sil.risk, y=CCND1))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=CaseID))+
  scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                "WU-07"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       y="CCND1 Z-score",
       color="Case", title = "CCND1")+
  theme(legend.position = "none", plot.title= element_text(face="bold", hjust=0.5))+
  stat_pvalue_manual(stat, y.position = c(3.5, 2.4, 3.1, 2.7), hide.ns = TRUE, label="p.adj.signif")

###Figure 2B
(p1+p2+p3)/(p4+p5+p6)

#ggsave("EPI6sharedgenes_0616.png", path=path)
p6+theme(legend.position="bottom")
#ggsave("legend_cases.png", path=path)


## Keratin 5 and Kertain 6
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

############Four gene risk score###################
risk.genes <- c("FN1", "KRT6A", "LAMC2", "TNC")
risk.data <- full.e[,colnames(full.e) %in% risk.genes]
risk.data$score <- rowSums(risk.data)/length(risk.genes)
full.e$score <- risk.data$score

stat <- wilcox_test(full.e, score~sil.risk)%>% add_significance()

ggplot(full.e, aes(x=sil.risk, y=score))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=CaseID))+
  scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                "WU-07"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="4 Gene Risk Score",
       title="4 Gene Risk Score")+
  stat_pvalue_manual(stat, y.position = c(3.1,1,2.8,2,3.3), label= "p.adj.signif", hide.ns = TRUE)+
  theme(plot.title=element_text(hjust=0.5, face="bold"),legend.position = "bottom")

####Figure 5A
#ggsave("4genescore_0616.png", path=path)

############### plots genes stroma ##################

full.s <- full[full$Compartment=="stroma",]

#stroma.genes ##id'ed by pathways
#"SFRP2" , "MMP9",   "OLR1"    "FN1" #CTHRC1 ###HK3 excluding case 18418
#"IL1RN"    "CHIT1" "HK3" ## not as good of a separation



stat <- wilcox_test(full.s, CD68~sil.risk)%>% add_significance()

p1 <- ggplot(full.s, aes(x=sil.risk, y=CD68))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=CaseID))+
  scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                "WU-07"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="CD68 Z-score",
       title="CD68")+
  theme(legend.position = "none", plot.title=element_text(hjust=0.5, face="bold"))+
  stat_pvalue_manual(stat, y.position=c(2.6,3.5,2.3,3,3.3,3.7), hide.ns=TRUE, lable="p.adj.signif")



stat <- wilcox_test(full.s, CTHRC1~sil.risk)%>% add_significance()

p2 <- ggplot(full.s, aes(x=sil.risk, y=CTHRC1))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=CaseID))+
  scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                "WU-07"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="CTHRC1 Z-score",
       title="CTHRC1")+
  theme(legend.position = "none", plot.title=element_text(hjust=0.5, face="bold"))+
  stat_pvalue_manual(stat, y.position = c(3.5,1,3.2,2.7,3.9), label= "p.adj.signif", hide.ns = TRUE)



stat <- wilcox_test(full.s, SFRP2~sil.risk)%>% add_significance()

p3 <- ggplot(full.s, aes(x=sil.risk, y=SFRP2))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=CaseID))+
  scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                "WU-07"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="SFRP2 Z-score",
       title="SFRP2")+
  theme(legend.position = "none", plot.title=element_text(hjust=0.5, face="bold"))+
  stat_pvalue_manual(stat, y.position = c(3.5,2.2,2.6,3), label= "p.adj.signif", hide.ns = TRUE)


stat <- wilcox_test(full.s, MMP9~sil.risk)%>% add_significance()

p4 <- ggplot(full.s, aes(x=sil.risk, y=MMP9))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=CaseID))+
  scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                "WU-07"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="MMP9 Z-score",
       title="MMP9")+ 
  theme(legend.position = "none", plot.title=element_text(hjust=0.5, face="bold"))+
  stat_pvalue_manual(stat, y.position = c(2.3,3.1,3.3), hide.ns = TRUE, label="p.adj.signif")

stat <- wilcox_test(full.s, FN1~sil.risk)%>% add_significance()

p5 <- ggplot(full.s, aes(x=sil.risk, y=FN1))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=CaseID))+
  scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                "WU-07"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="FN1 Z-score",
       title="FN1")+
  theme(legend.position = "none", plot.title=element_text(hjust=0.5, face="bold"))+
  stat_pvalue_manual(stat, y.position = c(2.3,3.5,1.7,2,2.7,2.8,3.1), hide.ns = TRUE, label="p.adj.signif")


stat <- wilcox_test(full.s, POSTN~sil.risk)%>% add_significance()

p6 <- ggplot(full.s, aes(x=sil.risk, y=POSTN))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=CaseID))+
  scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                "WU-07"= "#DEA0FD"))+
  scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                   labels=c("Benign", "A", "B Low", "B High", "C"))+
  labs(x=element_blank(),
       color="Case",
       y="POSTN Z-score",
       title="POSTN")+
  theme(legend.position = "none", plot.title=element_text(hjust=0.5, face="bold"))+
  stat_pvalue_manual(stat, y.position = c(3.6,3.2,3.4,3.8), hide.ns = TRUE, label="p.adj.signif")


####FIGURE 3B
(p1+p2+p3)/(p4+p5+p6)

#ggsave("selectStromaGenes2.png", path=path)






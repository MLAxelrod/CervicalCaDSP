###Cibersortx###

#run with LM22 with 100 permutations in absolute mode for job21
#run with LM22 with 1000 permutations in absolute mode for job22;l quantile normalization disable

library(plyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)

theme_set(theme_classic())
theme_update(axis.text=element_text(color="black", size=10), 
             axis.title.y = element_text(size=12), 
             axis.title.x = element_text(size=11))

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

data <-  read.delim("C:/Users/marga/Desktop/CCa_GEOMX/Cibersort/CIBERSORTx_Job22_Results.txt", sep="\t", row.names=1, check.names = TRUE)
rownames(data) <- paste0("X", rownames(data)) #R doesn't like row/col names starting with numbers

rownames(data)==rownames(meta)
check <- data.frame(x=rownames(data), y=rownames(meta), z=rownames(data)==rownames(meta)) ##all same, just formatted differently


full <- cbind(meta,data)
full.e <- full[full$Compartment=="epi",]
full.s <- full[full$Compartment=="stroma",]


imm <- colnames(full.s)[18:39]

df <- rbind(wilcox_test(full.s, B.cells.naive~sil.risk),
            wilcox_test(full.s, B.cells.memory~sil.risk),
            wilcox_test(full.s, Plasma.cells~sil.risk),
            wilcox_test(full.s, T.cells.CD8~sil.risk),
            wilcox_test(full.s, T.cells.CD4.naive~sil.risk),
            wilcox_test(full.s, T.cells.CD4.memory.resting~sil.risk),
            wilcox_test(full.s, T.cells.CD4.memory.activated~sil.risk),
            wilcox_test(full.s, T.cells.follicular.helper~sil.risk),
            wilcox_test(full.s, T.cells.regulatory..Tregs.~sil.risk),
            #wilcox_test(full.s, T.cells.gamma.delta~sil.risk), #too few to run
            wilcox_test(full.s, NK.cells.resting~sil.risk),
            wilcox_test(full.s, NK.cells.activated~sil.risk),
            wilcox_test(full.s, Monocytes~sil.risk),
            wilcox_test(full.s, Macrophages.M0~sil.risk),
            wilcox_test(full.s, Macrophages.M1~sil.risk),
            wilcox_test(full.s, Macrophages.M2~sil.risk),
            wilcox_test(full.s, Dendritic.cells.resting~sil.risk),
            wilcox_test(full.s, Dendritic.cells.activated~sil.risk),
            wilcox_test(full.s, Mast.cells.resting~sil.risk),
            wilcox_test(full.s, Mast.cells.activated~sil.risk),
            #wilcox_test(full.s, Eosinophils~sil.risk), #too few to run
            wilcox_test(full.s, Neutrophils~sil.risk))

df$BH <- p.adjust(df$p, method="BH")

#write.csv(df, file="C:/Users/marga/Desktop/CCa_GEOMX/Cibersort/stats2.csv")

immplt<-function(cells){
  ggplot(full.s, aes(x=sil.risk, y=cells))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(color=CaseID))+
    scale_color_manual(values = c("WU-01"="#F6222E", "WU-02"="#3283FE" , "WU-03"="#FEAF16", 
                                  "WU-04"="#C4451C","WU-05"="#2ED9FF", "WU-06"="#1C8356",
                                  "WU-07"= "#DEA0FD"))+
    scale_x_discrete(limits=c("normal.Normal", "A.Low", "B.Low", "B.High", "C.High"),
                     labels=c("Benign", "A", "B Low", "B High", "C"))+
    theme(legend.position = "none", plot.title= element_text(face="bold", hjust=0.5), 
          axis.title.y=element_text(size=14))
}





p1 <- immplt(full.s$B.cells.memory)+ylab("B cells memory")
p2 <- immplt(full.s$B.cells.naive)+ylab("B cells naive")
p3 <- immplt(full.s$Plasma.cells)+ylab("Plasma cells")
p4 <- immplt(full.s$T.cells.CD8)+ylab("T cells CD8")
p5 <- immplt(full.s$T.cells.CD4.naive)+ylab("T cells CD4 naive")
p6 <- immplt(full.s$T.cells.CD4.memory.resting)+ylab("T cells CD4 memory resting")
p7 <- immplt(full.s$T.cells.CD4.memory.activated)+ylab("T cells CD4 memory activated")
p8 <- immplt(full.s$T.cells.follicular.helper)+ylab("T cells follicular helper")
p9 <- immplt(full.s$T.cells.regulatory..Tregs.)+ylab("Tregs")
p10 <- immplt(full.s$T.cells.gamma.delta)+ylab("T cells gamma delta")
p11 <- immplt(full.s$NK.cells.resting)+ylab("NK cells resting")
p12 <- immplt(full.s$NK.cells.activated)+ylab("NK cells activated")
p13 <- immplt(full.s$Monocytes)+ylab("Monocytes")
p14 <- immplt(full.s$Macrophages.M0)+ylab("Macrophages M0")  
p15 <- immplt(full.s$Macrophages.M1)+ylab("Macrophages M1")  
p16 <- immplt(full.s$Macrophages.M2)+ylab("Macrophages M2")  
p17 <- immplt(full.s$Dendritic.cells.activated)+ylab("Dendritic cells activated")  
p18 <- immplt(full.s$Dendritic.cells.resting)+ylab("Dendritic cells resting")  
p19 <- immplt(full.s$Mast.cells.resting)+ylab("Mast cells resting")
p20 <- immplt(full.s$Mast.cells.activated)+ylab("Mast cells activated")  
p21 <- immplt(full.s$Eosinophils)+ylab("Eosinophils")
p22 <- immplt(full.s$Neutrophils)+ylab("Neutrophils")  

(p1|p2|p3|p4|p5)/(p6|p7|p8|p9|p10)/(p11|p12|p13|p14|p15|p16)/(p17|p18|p19|p20|p21|p22)
#ggsave("immplts.png", path=path)



full.s$B.cells <- full.s$B.cells.naive + full.s$B.cells.memory
full.s$T.cells <- full.s$T.cells.CD8 + full.s$T.cells.CD4.memory.activated + full.s$T.cells.CD4.memory.resting + full.s$T.cells.CD4.naive +full.s$T.cells.follicular.helper + full.s$T.cells.regulatory..Tregs.
#full.s$macs <- full.s$Macrophages.M0 + full.s$Macrophages.M1 + full.s$Macrophages.M2 + full.s$Monocytes
full.s$macs <- full.s$Macrophages.M0 + full.s$Macrophages.M1 + full.s$Macrophages.M2


stat <- wilcox_test(full.s, Monocytes~sil.risk)
p1 <- immplt(full.s$Monocytes)+
  labs(x=element_blank(),
       color="Case",
       y="Monocytes")+
  stat_pvalue_manual(stat, y.position = c(0.4), hide.ns = TRUE, label="p.adj.signif")


stat <- wilcox_test(full.s, Macrophages.M0~sil.risk)
p2 <- immplt(full.s$Macrophages.M0)+
  labs(x=element_blank(),
       color="Case",
       y="Macrophages.M0")+
  stat_pvalue_manual(stat, y.position = c(0.75), hide.ns = TRUE, label="p.adj.signif")


stat <- wilcox_test(full.s, Macrophages.M1~sil.risk)
p3 <- immplt(full.s$Macrophages.M1)+
  labs(x=element_blank(),
       color="Case",
       y="Macrophages.M1")+
  stat_pvalue_manual(stat, y.position = c(0.2), hide.ns = TRUE, label="p.adj.signif")

stat <- wilcox_test(full.s, Macrophages.M2~sil.risk)
p4 <- immplt(full.s$Macrophages.M1)+
  labs(x=element_blank(),
       color="Case",
       y="Macrophages.M2")+
  stat_pvalue_manual(stat, y.position = c(0.2, 0.215, 0.23), hide.ns = TRUE, label="p.adj.signif")


stat <- wilcox_test(full.s, T.cells~sil.risk)
p5 <- immplt(full.s$T.cells)+
  labs(x=element_blank(),
       color="Case",
       y="Total T cells")+
  stat_pvalue_manual(stat, y.position = c(2, 2.2), hide.ns = TRUE, label="p.adj.signif")


stat <- wilcox_test(full.s, B.cells~sil.risk)
p6 <- immplt(full.s$B.cells)+
  labs(x=element_blank(),
       color="Case",
       y="Total B cells")+
  stat_pvalue_manual(stat, y.position = c(0.5,0.55, 0.6), hide.ns = TRUE, label="p.adj.signif")


(p1+p2+p3)/(p4+p5+p6)

#ggsave("CibersortFig.png",path=path)


stat <- wilcox_test(full.s, macs~sil.risk)
immplt(full.s$macs)+
  labs(x=element_blank(),
       color="Case",
       y="Total Macrophages")+
  stat_pvalue_manual(stat, y.position = c(1.1, 1.15, 1.2), hide.ns = TRUE, label="p.adj.signif")

#ggsave("totalmacs.png", path=path)

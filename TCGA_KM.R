##TCGA KM 

library(survival)
library(survminer)
library(precrec)


scaleRow <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
  return(x)
}

path <- "C:/Users/marga/Desktop/CCa_GEOMX/plots/"

#####cbioportal data #####

data <- read_tsv("C:/Users/marga/Desktop/CCa_GEOMX/TCGA/cbio/mRNA expression (RNA Seq V2 RSEM).txt")
data <- read_tsv("C:/Users/marga/Desktop/CCa_GEOMX/TCGA/cbio/mRNA expression (RNA Seq V2 RSEM)_stroma.txt")
meta <- read_tsv("C:/Users/marga/Desktop/CCa_GEOMX/TCGA/cbio/KM_Plot__Overall_(months).txt")

colnames(meta) <- c("Study", "Pt", "OS_status", "OS_months")

##make status a numeric
meta$OS_stat_num <- as.numeric(substr(meta$OS_status,1,1))

##get rid of trailing "-01" in sample ID
data$SAMPLE_ID <- substr(data$SAMPLE_ID, 1, nchar(data$SAMPLE_ID)-3)

data$SAMPLE_ID %in% meta$Pt
all(data$SAMPLE_ID %in% meta$Pt)

df <- data[match(meta$Pt, data$SAMPLE_ID),]

all(df$SAMPLE_ID == meta$Pt)

##Scale Row by z-score
df2 <- as.data.frame(df[,3:177])
df2 <- as.data.frame(df[,3:15])
df2$SLURP2 <-  NULL
df2$TMEM265 <- NULL
sum(is.na(df2)) ##remove NA values
df4 <- as.data.frame(t(df2))

df3 <- scaleRow(df4)

df5 <- as.data.frame(t(df3))

summary(meta$OS_months)

meta <- mutate(meta, OSgroup=ifelse(meta$OS_months>22.57, "long", "short"))


full <- cbind(meta, df5)

cp <- cutpointr(full, KRT6A, OSgroup, method = maximize_metric, metric = youden, pos_class= "long", direction=">=")
summary(cp)
summary(full$KRT6A)


full <- mutate(full, KRT6Agroup=ifelse(full$KRT6A>-0.3751, "high", "low")) #youden
full <- mutate(full, KRT6Agroup=ifelse(full$KRT6A>-0.3428, "high", "low")) #median


fit<-survfit(Surv(OS_months, OS_stat_num) ~ KRT6Agroup, data=full)
print(fit)

ggsurvplot(fit,
           data=full,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)


###FN1
cp <- cutpointr(full, FN1, OSgroup, method = maximize_metric, metric = youden, pos_class= "long", direction=">=")
summary(cp)
summary(full$FN1)


full <- mutate(full, FN1group=ifelse(full$FN1>-0.3272, "high", "low")) #youden
full <- mutate(full, FN1group=ifelse(full$FN1>-0.4709, "high", "low")) #median
fit<-survfit(Surv(OS_months, OS_stat_num) ~ FN1group, data=full)

ggsurvplot(fit,
           data=full,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)


#LAMC2
cp <- cutpointr(full, LAMC2, OSgroup, method = maximize_metric, metric = youden, pos_class= "long")
summary(cp)
summary(full$LAMC2)
full <- mutate(full, LAMC2group=ifelse(full$LAMC2>-0.8430, "high", "low"))
full <- mutate(full, LAMC2group=ifelse(full$LAMC2>-0.3099, "high", "low"))

fit<-survfit(Surv(OS_months, OS_stat_num) ~ LAMC2group, data=full)
ggsurvplot(fit,
           data=full,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)
#p=0.17 for median

#TNC
cp <- cutpointr(full, TNC, OSgroup, method = maximize_metric, metric = youden, pos_class= "long")
summary(cp)
summary(full$TNC)
full <- mutate(full, TNCgroup=ifelse(full$TNC>-0.3100, "high", "low"))
full <- mutate(full, TNCgroup=ifelse(full$TNC>-0.7032, "high", "low"))
#p=0.25 for median ; p=0.19 for youden

fit<-survfit(Surv(OS_months, OS_stat_num) ~ TNCgroup, data=full)
ggsurvplot(fit,
           data=full,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)


####SERPINA3

cp <- cutpointr(full, SERPINA3, OSgroup, method = maximize_metric, metric = youden, pos_class= "long")
summary(cp)
summary(full$SERPINA3)
full <- mutate(full, SERPINA3group=ifelse(full$SERPINA3>-0.2429, "high", "low"))
full <- mutate(full, SERPINA3group=ifelse(full$SERPINA3>-0.2407, "high", "low"))
#p=0.25 for median ; p=0.19 for youden

fit<-survfit(Surv(OS_months, OS_stat_num) ~ SERPINA3group, data=full)
ggsurvplot(fit,
           data=full,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)



##CCND1
cp <- cutpointr(full, CCND1, OSgroup, method = maximize_metric, metric = youden, pos_class= "long")
summary(cp)
summary(full$CCND1)
full <- mutate(full, CCND1group=ifelse(full$CCND1>-0.5829, "high", "low"))
full <- mutate(full, CCND1group=ifelse(full$CCND1>-0.4275, "high", "low"))
#p=0.25 for median ; p=0.19 for youden

fit<-survfit(Surv(OS_months, OS_stat_num) ~ CCND1group, data=full)
ggsurvplot(fit,
           data=full,
           pval = TRUE, pval.coord=c(25,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)











####Make risk score 

#risk.genes <- c("SERPINA3", "CCDN1", "FN1", "KRT6A", "LAMC2", "TNC")
#risk.genes <- c("CCDN1", "FN1", "KRT6A", "LAMC2", "TNC")

risk.genes <- c("FN1", "KRT6A", "LAMC2", "TNC") ##epi genes

risk.genes <- c("FN1", "HK3", "MMP9", "CTHRC1", "SRFP2", "OLR1") ##stroma

risk.data <- full[,colnames(full) %in% risk.genes]

risk.data$score <- rowSums(risk.data)/length(risk.genes)

full$score <- risk.data$score


cp <- cutpointr(full, score, OSgroup, method = maximize_metric, metric = youden, pos_class= "long")
summary(cp)
summary(full$score)

#six gene
#full <- mutate(full, SCOREgroup=ifelse(full$score>-0.0019, "high", "low"))
#full <- mutate(full, SCOREgroup=ifelse(full$score>-0.09583, "high", "low"))


#four gene

full <- mutate(full, SCOREgroup=ifelse(full$score>-0.3076, "high", "low"))
full <- mutate(full, SCOREgroup=ifelse(full$score>-0.1380, "high", "low"))


#stroma 
full <- mutate(full, SCOREgroup=ifelse(full$score>0.3189, "high", "low"))
full <- mutate(full, SCOREgroup=ifelse(full$score>-0.1200, "high", "low"))



fit<-survfit(Surv(OS_months, OS_stat_num) ~ SCOREgroup, data=full)
ggsurvplot(fit,
           data=full,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)



foo <- cbind(full$OS_months, full$OS_status, full$score, full$SCOREgroup)
#for sanity check that association all in expected directions

#######






################ Stroma genes
#"SFRP2" , "MMP9",   "OLR1"    "FN1" #CTHRC1 ###HK3 excluding case 18418; already did FN1

###MMP9 ## No association
cp <- cutpointr(full, MMP9, OSgroup, method = maximize_metric, metric = youden, pos_class= "long", direction=">=")
summary(cp)
summary(full$MMP9)
full <- mutate(full, Genegroup=ifelse(full$MMP9>0.442, "high", "low")) #youden
full <- mutate(full, Genegroup=ifelse(full$MMP9>-0.2337, "high", "low")) #median

#SFRP9 #no association
cp <- cutpointr(full, SFRP2, OSgroup, method = maximize_metric, metric = youden, pos_class= "long", direction=">=")
summary(cp)
summary(full$SFRP2)
full <- mutate(full, Genegroup=ifelse(full$SFRP2>0.0554, "high", "low")) #youden
full <- mutate(full, Genegroup=ifelse(full$SFRP2>-0.2337, "high", "low")) #median


#HK3 #no association
cp <- cutpointr(full, HK3, OSgroup, method = maximize_metric, metric = youden, pos_class= "long", direction=">=")
summary(cp)
summary(full$HK3)
full <- mutate(full, Genegroup=ifelse(full$HK3>-0.1166, "high", "low")) #youden
full <- mutate(full, Genegroup=ifelse(full$HK3>-0.32164, "high", "low")) #median


#CTHRC1 #no assoicaiton
cp <- cutpointr(full, CTHRC1, OSgroup, method = maximize_metric, metric = youden, pos_class= "long", direction=">=")
summary(cp)
summary(full$CTHRC1)
full <- mutate(full, Genegroup=ifelse(full$CTHRC1>0.2916, "high", "low")) #youden
full <- mutate(full, Genegroup=ifelse(full$CTHRC1>-0.3976, "high", "low")) #median


#OLR1 # no association
cp <- cutpointr(full, OLR1, OSgroup, method = maximize_metric, metric = youden, pos_class= "long", direction=">=")
summary(cp)
summary(full$OLR1)
full <- mutate(full, Genegroup=ifelse(full$OLR1>-0.0428, "high", "low")) #youden
full <- mutate(full, Genegroup=ifelse(full$OLR1>-0.4577, "high", "low")) #median




fit<-survfit(Surv(OS_months, OS_stat_num) ~ Genegroup, data=full)
ggsurvplot(fit,
           data=full,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)















meta$SERPINA3 <- df$SERPINA3

cp <- cutpointr(meta, SERPINA3, OSgroup, method = maximize_metric, metric = youden, pos_class= "long", direction=">=")
summary(cp)
summary(meta$SERPINA3)

meta <- mutate(meta, SERPgroup=ifelse(meta$SERPINA3>602, "high", "low"))
meta <- mutate(meta, SERPgroup=ifelse(meta$KRT6A>642, "high", "low"))


fit<-survfit(Surv(OS_months, OS_stat_num) ~ SERPgroup, data=meta)
print(fit)

ggsurvplot(fit,
           data=meta,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)


meta$CCND1 <- df$CCND1
cp <- cutpointr(meta, CCND1, OSgroup, method = maximize_metric, metric = youden, pos_class= "long", direction=">=")
summary(cp)
summary(meta$CCND1)
meta <- mutate(meta, GENEgroup=ifelse(meta$CCND1>527, "high", "low"))


fit<-survfit(Surv(OS_months, OS_stat_num) ~ GENEgroup, data=meta)

ggsurvplot(fit,
           data=meta,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)

##FN1
meta$FN1 <- df$FN1
summary(meta$FN1)
meta <- mutate(meta, GENEgroup=ifelse(meta$FN1>2025, "high", "low"))


fit<-survfit(Surv(OS_months, OS_stat_num) ~ GENEgroup, data=meta)

ggsurvplot(fit,
           data=meta,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)

##p = 0.26 , most promising so far ....

















###tcga data ### but questionable normalization
data <- read_tsv("C:/Users/marga/Desktop/CCa_GEOMX/TCGA/HierCluster.2024-12-09.tsv")
meta <- read_tsv("C:/Users/marga/Desktop/CCa_GEOMX/TCGA/overall-survival-plot.2024-12-09.tsv")

genes <- colnames(data)

"POSTN" %in% genes

colnames(meta) <- c("id", "days", "months", "years", "censored", "survivalEstimate", "submitter_id", "project_id")

fit<-survfit(Surv(days, censored) ~ project_id, data=meta)
print(fit)

ggsurvplot(fit,
           data=meta,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)

summary(data$KRT6A)


data.filt <- data[data$Case %in% meta$submitter_id,]

df <- data.filt[match(meta$submitter_id, data.filt$Case),]

all(df$Case==meta$submitter_id)



meta$KRT6A <- df$KRT6A

summary(meta$days)

meta <- mutate(meta, OSgroup=ifelse(meta$days>505, "long", "short"))

cp <- cutpointr(meta, KRT6A, OSgroup, method = maximize_metric, metric = youden, pos_class= "long", direction=">=")
summary(cp)
summary(meta$KRT6A)

meta <- mutate(meta, KRT6Agroup=ifelse(meta$KRT6A>-0.32083, "high", "low"))


fit<-survfit(Surv(days, censored) ~ KRT6Agroup, data=meta)
print(fit)

ggsurvplot(fit,
           data=meta,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)


### not sure how the data you downloaded are normalized, only 34 have survival data, try cBioportal and see if you get the same results 




meta$POSTN <- df$POSTN


cp <- cutpointr(meta, POSTN, OSgroup, method = maximize_metric, metric = youden, pos_class= "long", direction=">=")
summary(cp)
summary(meta$POSTN)

meta <- mutate(meta, POSTNgroup=ifelse(meta$POSTN>-0.23948, "high", "low"))


fit<-survfit(Surv(days, censored) ~ POSTNgroup, data=meta)
print(fit)

ggsurvplot(fit,
           data=meta,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE)


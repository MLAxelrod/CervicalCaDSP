##TCGA KM 

library(plyr)
library(broom)
library(tidyverse)
library(survival)
library(survminer)
library(precrec)
library(cutpointr)
library(ggfortify)


scaleRow <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
  return(x)
}

path <- "C:/Users/marga/Desktop/CCa_GEOMX/plots/"

#####cbioportal data #####
### all data downloaded from cbioportal with differing genes included

data <- read_tsv("C:/Users/marga/Desktop/CCa_GEOMX/TCGA/cbio/mRNA expression (RNA Seq V2 RSEM).txt")
#data <- read_tsv("C:/Users/marga/Desktop/CCa_GEOMX/TCGA/cbio/mRNA expression (RNA Seq V2 RSEM)_stroma.txt")
#data <- read_tsv("C:/Users/marga/Desktop/CCa_GEOMX/TCGA/cbio/mRNA expression (RNA Seq V2 RSEM)_CD68.txt")
meta <- read_tsv("C:/Users/marga/Desktop/CCa_GEOMX/TCGA/cbio/KM_Plot__Overall_(months).txt")
stage <- read_tsv("C:/Users/marga/Desktop/CCa_GEOMX/TCGA/cbio/American_Joint_Committee_on_Cancer_Tumor_Stage_Code.full.txt")

colnames(meta) <- c("Study", "Pt", "OS_status", "OS_months")
colnames(stage) <- c("Study", "Pt", "AJCC_stage")


##make status a numeric
meta$OS_stat_num <- as.numeric(substr(meta$OS_status,1,1))

##add stage to meta
stage$Pt %in% meta$Pt
stage.df <- stage[match(meta$Pt, stage$Pt),]
stage.df$AJCC_stage2 <- replace_na(stage.df$AJCC_stage, "Unknown")
stage.df$AJCC_stage2 <- revalue(stage.df$AJCC_stage2, c("TX" = "Unknown"))
stage.df$Pt == meta$Pt
meta$Stage <- stage.df$AJCC_stage2
table(meta$Stage)


##get rid of trailing "-01" in sample ID
data$SAMPLE_ID <- substr(data$SAMPLE_ID, 1, nchar(data$SAMPLE_ID)-3)

##Match meta to data##MAJCC_stageatch meta to data
data$SAMPLE_ID %in% meta$Pt
all(data$SAMPLE_ID %in% meta$Pt)

df <- data[match(meta$Pt, data$SAMPLE_ID),]

all(df$SAMPLE_ID == meta$Pt)



##Scale Row by z-score; depending on number of genes included
df2 <- as.data.frame(df[,3:177])
#df2 <- as.data.frame(df[,3:15])
#df2 <- as.data.frame(df[,3:3])

#remove NAs
df2$SLURP2 <-  NULL
df2$TMEM265 <- NULL
sum(is.na(df2)) 

df4 <- as.data.frame(t(df2))
df3 <- scaleRow(df4)
df5 <- as.data.frame(t(df3))

summary(meta$OS_months)

meta <- mutate(meta, OSgroup=ifelse(meta$OS_months>22.57, "long", "short"))


full <- cbind(meta, df5)


####Make risk score 

risk.genes <- c("FN1", "KRT6A", "LAMC2", "TNC") ##epi genes

risk.data <- full[,colnames(full) %in% risk.genes]

risk.data$score <- rowSums(risk.data)/length(risk.genes)

full$score <- risk.data$score


cp <- cutpointr(full, score, OSgroup, method = maximize_metric, metric = youden, pos_class= "long")
summary(cp)
summary(full$score)


#four gene

full <- mutate(full, SCOREgroup=ifelse(full$score>-0.3076, "high", "low")) #by youden
#full <- mutate(full, SCOREgroup=ifelse(full$score>-0.1380, "high", "low")) #by median




fit<-survfit(Surv(OS_months, OS_stat_num) ~ SCOREgroup, data=full)
ggsurvplot(fit,
           data=full,
           pval = TRUE, pval.coord=c(20,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE,
           tables.col="black",
           palette = c("#DC3220", "#005AB5"))



table(full$SCOREgroup, full$Stage)

full$stage_simple <- revalue(full$Stage, c("T1b" = "T1", "T1b1" = "T1", "T1b2" = "T1", "T2a" = "T2", "T2a1" = "T2", "T2a2" = "T2", "T2b" = "T2", "T3b" = "T3"))
table(full$stage_simple, full$SCOREgroup)

fit<-survfit(Surv(OS_months, OS_stat_num) ~ stage_simple, data=full)
ggsurvplot(fit,
           data=full,
           pval = TRUE, pval.coord=c(20,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE,
           tables.col="black")


##############################
### COX proportional hazard model

#cox <- coxph(Surv(OS_months, OS_stat_num) ~ score +strata(stage_simple), data=full)

cox <- coxph(Surv(OS_months, OS_stat_num) ~ score, data=full)

summary(cox)

tcox <- tidy(cox, exponentiate=T, conf.int=T) 

cox_fit <- survfit(cox)
#plot(cox_fit, main = "cph model", xlab="Time")
autoplot(cox_fit)



#########################################################################################################
####Additional single gene analyses#########

###CD68
cp <- cutpointr(full, CD68, OSgroup, method = maximize_metric, metric = youden, pos_class= "long", direction=">=")
summary(cp)
summary(full$CD68)


full <- mutate(full, CD68group=ifelse(full$CD68>-0.921, "high", "low")) #youden
full <- mutate(full, CD68group=ifelse(full$CD68>-0.2353, "high", "low")) #median


fit<-survfit(Surv(OS_months, OS_stat_num) ~ CD68group, data=full)
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




########## KRT6A
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




#############################
foo <- cbind(full$OS_months, full$OS_status, full$score, full$SCOREgroup)
#for sanity check that association all in expected directions


################ Stroma genes
#"SFRP2" , "MMP9",   "OLR1"    "FN1" #CTHRC1 ###HK3 

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



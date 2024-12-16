##TCGA KM 

library(survival)
library(survminer)
library(precrec)

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


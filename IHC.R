
ihc <- read.csv("C:/Users/marga/Desktop/CCa_GEOMX/data/IHCdata.csv")

path <- "C:/Users/marga/Desktop/CCa_GEOMX/plots/"

table(ihc$SILVA)

sub <- ihc[ihc$SILVA %in% c("A", "B", "C", "Benign"),]

sub2 <- sub[!sub$CD68.MA=="n",]
sub2 <- sub2[!sub2$CK56.I.MA=="n",]

sub2$CK56.P.LS <- as.numeric(sub2$CK56.P.LS)
sub2$CK56.P.MA <- as.numeric(sub2$CK56.P.MA)

sub2$CK56.I.MA <- as.numeric(sub2$CK56.I.MA)
sub2$CK56.I.LS <- as.numeric(sub2$CK56.I.LS)

sub2$CD68.LS <- as.numeric(sub2$CD68.LS)
sub2$CD68.MA <- as.numeric(sub2$CD68.MA)

p1 <- ggplot(sub2, aes(CK56.I.LS, CK56.I.MA))+
  geom_jitter(width=0.1, height=0.1)+
  labs(x="CK5/6 intensity LS",
       y="CK5/6 intensity MA")



p2 <- ggplot(sub2, aes(CK56.P.LS, CK56.P.MA))+
  geom_point()+
  labs(x="CK5/6 percent LS",
       y="CK5/6 percent MA")
  

p3 <- ggplot(sub2, aes(CD68.LS, CD68.MA))+
  geom_jitter(width=0.1, height=0.1)+
  labs(x="CD68 score LS", 
       y="CD68 score MA")



p1+p2+p3
#ggsave("concordance.png", path=path)



###calculate averages

sub2$avgCD68 <- (sub2$CD68.MA+sub2$CD68.LS)/2
sub2$avgCK56.I <- (sub2$CK56.I.MA+sub2$CK56.I.LS)/2
sub2$avgCK56.P <- (sub2$CK56.P.MA+sub2$CK56.P.LS)/2

sub2$CK56.C <- (sub2$avgCK56.I*sub2$avgCK56.P)


stat <- wilcox_test(sub2, avgCD68~SILVA)%>% add_significance()

ggplot(sub2, aes(SILVA, avgCD68))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1, height=0)+
  scale_x_discrete(limits=c("Benign", "A", "B", "C"))+
  stat_pvalue_manual(stat, y.position=c(1.5, 2.4, 2.2, 2.5), hide.ns=TRUE, lable="p.adj.signif")+
  labs(x="Silva Pattern",
       y="CD68 average score",
         title="CD68 IHC")+
  theme(plot.title = element_text(hjust=0.5))

#ggsave("CD68avgscore_s.png", path=path)


stat <- wilcox_test(sub2, avgCK56.I~SILVA)%>% add_significance()

p1 <- ggplot(sub2, aes(SILVA, avgCK56.I))+
 geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1, height=0)+
  scale_x_discrete(limits=c("Benign", "A", "B", "C"))+
  stat_pvalue_manual(stat, y.position=c(2.3), hide.ns=TRUE, lable="p.adj.signif", step.increase = 0.1)+
  labs(x="Silva Pattern",
       y="CK5/6 intensity average score",
       title="CK5/6 IHC intensity")+
  theme(plot.title = element_text(hjust=0.5))



stat <- wilcox_test(sub2, avgCK56.P~SILVA)%>% add_significance()

p2 <- ggplot(sub2, aes(SILVA, avgCK56.P))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1, height=0)+
  scale_x_discrete(limits=c("Benign", "A", "B", "C"))+
  stat_pvalue_manual(stat, y.position=c(101), hide.ns=TRUE, lable="p.adj.signif", step.increase = 0.1)+
  labs(x="Silva Pattern",
       y="CK5/6 percent average",
       title="CK5/6 IHC percent")+
  theme(plot.title = element_text(hjust=0.5))


stat <- wilcox_test(sub2, CK56.C~SILVA)%>% add_significance()

p3 <- ggplot(sub2, aes(SILVA, CK56.C))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1, height=0)+
  scale_x_discrete(limits=c("Benign", "A", "B", "C"))+
  stat_pvalue_manual(stat, y.position=c(201), hide.ns=TRUE, lable="p.adj.signif", step.increase = 0.1)+
  labs(x="Silva Pattern",
       y="CK5/6 (percent * intensity)",
       title="CK5/6 IHC composite")+
  theme(plot.title = element_text(hjust=0.5))

p3
#ggsave("ck56_c2.png", path=path)

p1+p2+p3
#ggsave("CK56ihc.png", path=path)


####
stat <- wilcox_test(sub2, CD68.LS~SILVA)%>% add_significance()
p1 <- ggplot(sub2, aes(SILVA, CD68.LS))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1, height=0)+
  scale_x_discrete(limits=c("Benign", "A", "B", "C"))+
  stat_pvalue_manual(stat, y.position=c(2), hide.ns=TRUE, lable="p.adj.signif", step.increase = 0.1)+
  labs(x="Silva Pattern",
       y="CD68 LS score",
       title="CD68 IHC")+
  theme(plot.title = element_text(hjust=0.5))


stat <- wilcox_test(sub2, CD68.MA~SILVA)%>% add_significance()
p2 <- ggplot(sub2, aes(SILVA, CD68.MA))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1, height=0)+
  scale_x_discrete(limits=c("Benign", "A", "B", "C"))+
  stat_pvalue_manual(stat, y.position=c(2.2), hide.ns=TRUE, lable="p.adj.signif", step.increase = 0.1)+
  labs(x="Silva Pattern",
       y="CD68 MA score",
       title="CD68 IHC")+
  theme(plot.title = element_text(hjust=0.5))

p1+p2
#ggsave("CD68_each.png", path=path)

c <- sub2[sub2$SILVA=="C",]
table(sub2$Case)

sub3 <- sub2[!sub2$Case=="18418",]


stat <- wilcox_test(sub3, avgCD68~SILVA)%>% add_significance()

p1 <- ggplot(sub3, aes(SILVA, avgCD68))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1, height=0)+
  scale_x_discrete(limits=c("Benign", "A", "B", "C"))+
  stat_pvalue_manual(stat, y.position=c(1.5, 2.4, 2.2, 2.55, 2.65), hide.ns=TRUE, lable="p.adj.signif")+
  labs(x="Silva Pattern",
       y="CD68 average score",
       title="CD68 IHC")+
  theme(plot.title = element_text(hjust=0.5))

stat <- wilcox_test(sub3, CK56.C~SILVA)%>% add_significance()

p2 <- ggplot(sub3, aes(SILVA, CK56.C))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1, height=0)+
  scale_x_discrete(limits=c("Benign", "A", "B", "C"))+
  stat_pvalue_manual(stat, y.position=c(201), hide.ns=TRUE, lable="p.adj.signif", step.increase = 0.1)+
  labs(x="Silva Pattern",
       y="CK5/6 (percent * intensity)",
       title="CK5/6 IHC composite")+
  theme(plot.title = element_text(hjust=0.5))

(p1+p2)+plot_annotation(title="Without case 18418", theme=theme(plot.title = element_text(hjust=0.5, size=15)))
ggsave("ihc_wo_18.png", path=path)

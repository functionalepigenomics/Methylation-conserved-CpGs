setwd("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/Mdata/Genotype-Methylation/Figure/Beta_correlation")
MC_CpG<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/MSCC.txt")
SC_CpG<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/SCC.txt")
NC_CpG<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/NCC.txt")

zz=gzfile('~/KoborLab/kobor_space/zdong/Monkey/Probes3/Version2-Genotype-Methylation/Brain_GSE112525/Correlation_MAF0.05/FDR/MHB-SNP_FDRlessthan0.05.txt.gz','rt')
re<-read.table(zz,header = T)
re1<-re
brain_mc<-abs(re1$beta[re1$gene%in%MC_CpG$V1])
brain_mc<-data.frame(brain_mc,"MCC")
names(brain_mc)<-c("mpg","cyl")
brain_sc<-abs(re1$beta[re1$gene%in%SC_CpG$V1])
brain_sc<-data.frame(brain_sc,"SCC")
names(brain_sc)<-c("mpg","cyl")
brain_nc<-abs(re1$beta[re1$gene%in%NC_CpG$V1])
brain_nc<-data.frame(brain_nc,"NCC")
names(brain_nc)<-c("mpg","cyl")
brain<-rbind(brain_mc,brain_sc,brain_nc)
brain$tissue<-"Brain"

zz=gzfile('~/KoborLab/kobor_space/zdong/Monkey/Probes3/Version2-Genotype-Methylation/Whole-blood/Correlation_MAF0.05/FDR/MHB-SNP_FDRlessthan0.05.txt.gz','rt')
re<-read.table(zz,header = T)
re1<-re
blood_mc<-abs(re1$beta[re1$gene%in%MC_CpG$V1])
blood_mc<-data.frame(blood_mc,"MCC")
names(blood_mc)<-c("mpg","cyl")
blood_sc<-abs(re1$beta[re1$gene%in%SC_CpG$V1])
blood_sc<-data.frame(blood_sc,"SCC")
names(blood_sc)<-c("mpg","cyl")
blood_nc<-abs(re1$beta[re1$gene%in%NC_CpG$V1])
blood_nc<-data.frame(blood_nc,"NCC")
names(blood_nc)<-c("mpg","cyl")
blood<-rbind(blood_mc,blood_sc,blood_nc)
blood$tissue<-"Blood"

zz=gzfile('~/KoborLab/kobor_space/zdong/Monkey/Probes3/Version2-Genotype-Methylation/Buccal_GECKO/Correlation_MAF0.05/Lastest/FDR/MHB-SNP_FDRlessthan0.05.txt.gz','rt')
re<-read.table(zz,header = T)
re1<-re
buccal_mc<-abs(re1$beta[re1$gene%in%MC_CpG$V1])
buccal_mc<-data.frame(buccal_mc,"MCC")
names(buccal_mc)<-c("mpg","cyl")
buccal_sc<-abs(re1$beta[re1$gene%in%SC_CpG$V1])
buccal_sc<-data.frame(buccal_sc,"SCC")
names(buccal_sc)<-c("mpg","cyl")
buccal_nc<-abs(re1$beta[re1$gene%in%NC_CpG$V1])
buccal_nc<-data.frame(buccal_nc,"NCC")
names(buccal_nc)<-c("mpg","cyl")
buccal<-rbind(buccal_mc,buccal_sc,buccal_nc)
buccal$tissue<-"Buccal"

zz=gzfile('~/KoborLab/kobor_space/zdong/Monkey/Probes3/Version2-Genotype-Methylation/LCLs_GSE24277//Correlation_MAF0.05/FDR/MHB-SNP_FDRlessthan0.05.txt.gz','rt')
re<-read.table(zz,header = T)
re1<-re
LCL_mc<-abs(re1$beta[re1$gene%in%MC_CpG$V1])
LCL_mc<-data.frame(LCL_mc,"MCC")
names(LCL_mc)<-c("mpg","cyl")
LCL_sc<-abs(re1$beta[re1$gene%in%SC_CpG$V1])
LCL_sc<-data.frame(LCL_sc,"SCC")
names(LCL_sc)<-c("mpg","cyl")
LCL_nc<-abs(re1$beta[re1$gene%in%NC_CpG$V1])
LCL_nc<-data.frame(LCL_nc,"NCC")
names(LCL_nc)<-c("mpg","cyl")
LCL<-rbind(LCL_mc,LCL_sc,LCL_nc)
LCL$tissue<-"LCLs"


zz=gzfile('~/KoborLab/kobor_space/zdong/Monkey/Probes3/Version2-Genotype-Methylation/Saliva_GSE99091/Correlation_MAF0.05/FDR/MHB-SNP_FDRlessthan0.05.txt.gz','rt')
re<-read.table(zz,header = T)
re1<-re
Saliva_mc<-abs(re1$beta[re1$gene%in%MC_CpG$V1])
Saliva_mc<-data.frame(Saliva_mc,"MCC")
names(Saliva_mc)<-c("mpg","cyl")
Saliva_sc<-abs(re1$beta[re1$gene%in%SC_CpG$V1])
Saliva_sc<-data.frame(Saliva_sc,"SCC")
names(Saliva_sc)<-c("mpg","cyl")
Saliva_nc<-abs(re1$beta[re1$gene%in%NC_CpG$V1])
Saliva_nc<-data.frame(Saliva_nc,"NCC")
names(Saliva_nc)<-c("mpg","cyl")
Saliva<-rbind(Saliva_mc,Saliva_sc,Saliva_nc)
Saliva$tissue<-"Saliva"

zz=gzfile('~/KoborLab/kobor_space/zdong/Monkey/Probes3/Version2-Genotype-Methylation/Skin_GSE53261/Correlation_MAF0.05/FDR/MHB-SNP_FDRlessthan0.05.txt.gz','rt')
re<-read.table(zz,header = T)
re1<-re
Skin_mc<-abs(re1$beta[re1$gene%in%MC_CpG$V1])
Skin_mc<-data.frame(Skin_mc,"MCC")
names(Skin_mc)<-c("mpg","cyl")
Skin_sc<-abs(re1$beta[re1$gene%in%SC_CpG$V1])
Skin_sc<-data.frame(Skin_sc,"SCC")
names(Skin_sc)<-c("mpg","cyl")
Skin_nc<-abs(re1$beta[re1$gene%in%NC_CpG$V1])
Skin_nc<-data.frame(Skin_nc,"NCC")
names(Skin_nc)<-c("mpg","cyl")
Skin<-rbind(Skin_mc,Skin_sc,Skin_nc)
Skin$tissue<-"Skin"
require(ggpubr)
library(RColorBrewer)
dfm<-rbind(blood,brain,buccal,LCL,Saliva,Skin)
dfm$cyl <- factor(dfm$cyl, levels = c("MCC", "SCC", "NCC"))
theme_set(theme_bw())
p<-ggplot(dfm, aes(x=cyl, y=mpg,fill = factor(cyl))) + 
  geom_boxplot(outlier.shape = NA) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=14), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text = element_text(size=14),
        legend.title =element_blank()
  )+ylab("Average regression coefficient for SNP-CpG pairs")+scale_y_continuous(limits = c(0, 0.4))+
  scale_fill_manual(values=brewer.pal(3,"Dark2"))+
  facet_wrap(~tissue)
ggsave("beta-allSNP_NEW.pdf",p, width=3, height=2, units="in", scale=3)

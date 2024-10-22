setwd("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/Mdata/Monkey_probes_selection/PhastCons-scores")
library(ggplot2)
require(ggridges)
library(dplyr)
library(RColorBrewer)
MC<-read.table('~/KoborLab/kobor_space/zdong/Monkey/Probes3/PhastCons-scores/MC-CpG.phastCons.txt.gz',sep='\t')
SC<-read.table('~/KoborLab/kobor_space/zdong/Monkey/Probes3/PhastCons-scores/SC-CpG.phastCons.txt.gz',sep='\t')
NC<-read.table('~/KoborLab/kobor_space/zdong/Monkey/Probes3/PhastCons-scores/NC-CpG.phastCons_removeXandY.txt.gz',sep='\t')
data<-rbind(MC,SC,NC)
data$V1<-gsub("chr", "", data$V1)
data$V2<-as.numeric(data$V2)+1
data$V3<-as.numeric(data$V2)+1
rownames(data)<-paste(data$V1,data$V2,data$V3,sep=":")
MSCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/MSCC_pos.bed")
SCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/SCC_pos.bed")
NCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/NCC_pos.bed")
rownames(MSCC)<-paste(MSCC$V1,MSCC$V2,MSCC$V3,sep=":")
rownames(SCC)<-paste(SCC$V1,SCC$V2,SCC$V3,sep=":")
rownames(NCC)<-paste(NCC$V1,NCC$V2,NCC$V3,sep=":")
MC<-data[rownames(data) %in% rownames(MSCC),]
SC<-data[rownames(data) %in% rownames(SCC),]
NC<-data[rownames(data) %in% rownames(NCC),]
MC$name<-1
SC$name<-0
NC$name<-2
MC<-MC[,5:6]
SC<-SC[,5:6]
NC<-NC[,5:6]
dat<-rbind(MC,SC,NC)
dat$V5<-as.numeric(dat$V5)
mean(dat$V5[dat$name==1],na.rm=TRUE) #0.361413
mean(dat$V5[dat$name==0],na.rm=TRUE) #0.2551309
mean(dat$V5[dat$name==2],na.rm=TRUE) #0.1245935
dat<-na.omit(dat)
# 1000 permutations
n=0
for (i in 1:1000){
  random<-sample(dat$V5[dat$name==0], size=length(dat$V5[dat$name==1]), replace=F)
  if (mean(random)>=mean(dat$V5[dat$name==1])){n=n+1}
}
n/1000 # 0

n=0
for (i in 1:1000){
  random<-sample(dat$V5[dat$name==2], size=length(dat$V5[dat$name==1]), replace=F)
  if (mean(random)>=mean(dat$V5[dat$name==1])){n=n+1}
}
n/1000 # 0

n=0
for (i in 1:1000){
  random<-sample(dat$V5[dat$name==2], size=length(dat$V5[dat$name==0]), replace=F)
  if (mean(random)>=mean(dat$V5[dat$name==0])){n=n+1}
}
n/1000 # 0

## ====== boxplot
t1<-dat$V5[dat$name==1]
t1<-data.frame(t1,"MCC")#"Population-specific effect and allelic frequency")
colnames(t1)<-c("Pst",'Type')
t2<-dat$V5[dat$name==0]
t2<-data.frame(t2,'SCC')#"Population-specific effect not allelic frequency")
colnames(t2)<-c("Pst",'Type')
t3<-dat$V5[dat$name==2]
t3<-data.frame(t3,'NCC')#"Population-specific effect not allelic frequency")
colnames(t3)<-c("Pst",'Type')
data<-rbind(t3,t2,t1)
colnames(data)<-c("Pst",'Type')
data$Type <- factor(data$Type, levels = c("MCC", "SCC","NCC"))

library(readr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)
raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  #axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
g <- ggplot(data, aes(y = Pst, x = Type, fill =Type)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  # geom_point(aes(y = Pst, color = Type),position = position_jitter(width = .15), size = .5, alpha = 0.1) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.8) +
  stat_summary(fun.y = mean, geom = "point", shape=20, size=3, color="black", fill="black",
               position = position_dodge(1.1))+
  expand_limits(x = 3) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  #scale_color_brewer(palette = "Set1") +
  #scale_fill_brewer(palette = "Set1") +
  scale_fill_manual(values=brewer.pal(3,"Dark2"))+
  scale_color_manual(values=brewer.pal(3,"Dark2"))+
  # coord_flip() +
  theme_bw() +
  raincloud_theme+theme(panel.grid = element_blank(),legend.position = "none", 
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        # panel.background = element_rect(fill = 'transparent', color = 'black'), 
                        legend.title = element_blank(), legend.key = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        strip.text = element_text(size=14), 
                        axis.title.x=element_blank(),
                        axis.title.y=element_text(size=14),
                        axis.text = element_text(size=14))+ ylab("PhastCons score")
ggsave("PhastCons_comparsion.pdf",g, width=5, height=3)#, units="in", scale=3)

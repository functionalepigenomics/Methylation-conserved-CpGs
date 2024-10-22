setwd("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/Mdata/Human_population/Tissue_figure/Boxplot")
library(ggpubr)		
library(RColorBrewer)	
library(ggplot2)
library(stringr)
theme_set(theme_classic())		
mpg<-read.table("../total1.log.gz",header=F,sep="\t")		
# # remove Altitude	
# mpg1<-mpg[mpg$V4 != "Altitude",]	
# mpg<-mpg1	

colnames(mpg)<-c("ID","Methylation Variance","con","Cancer")		
mpg$`Methylation Variance`<-abs(mpg$`Methylation Variance`)		
yy<-mpg$`Methylation Variance`		
xx<-mpg$Cancer	
zz<-mpg$con	
mpg<-data.frame(xx,yy,zz)	
mpg$zz[mpg$zz == 1]<-'MCC'	
mpg$zz[mpg$zz == 0]<-'SCC'	
mpg$zz[mpg$zz == 2]<-'NCC'	
# Plot		
## chang the order of x-axis		
mpg$xx <- factor(mpg$xx, levels = unique(mpg$xx))	#levels = rev(unique(mpg$xx)))		
mpg$zz <- factor(mpg$zz, levels = c('MCC','SCC','NCC'))	#rev(
mpg<-mpg[!is.na(mpg$zz),]	
mpg<-mpg[!is.na(mpg$yy),]	
# source('summarySE.R')	
#tgc<-summarySE(mpg, measurevar="yy", groupvars=c("zz","xx"))	
mpg1<-mpg[mpg$xx %in% c("Blood","Buccal","Cord-blood","Placenta","Sperm"),]
# mpg1 <- data.frame(lapply(mpg1, function(x) {gsub("Cord-blood", "Cord blood", x)}))

g <- ggplot(mpg1, aes(x=zz, y=yy, fill=zz)) + 
  # geom_boxplot(aes(middle = mean(yy))) +
  geom_boxplot(fatten = 2,outlier.shape = NA,position=position_dodge(1.1))+
  stat_summary(fun.y = mean, geom = "point", shape=20, size=3, color="black", fill="black",
               position = position_dodge(1.1))+
  facet_wrap(~xx, scale="free",nrow=2,labeller = as_labeller(c(
    "Blood"="Blood","Buccal"="Buccal",
    "Cord-blood"="Cord blood",
    "Placenta"="Placenta","Sperm"="Sperm")))+
  scale_fill_manual(values = brewer.pal(3,"Dark2"))+
  xlab("") +	
  ylab("Methylation variability")+ coord_cartesian(ylim = c(0, 0.2))+
  labs(fill="CpG type")+theme(legend.position = "none",
                              axis.title = element_text(size = 12),	strip.text = element_text(size = 14),
                              axis.text = element_text(size = 12),panel.spacing.x = unit(2, "lines"),
                              panel.spacing.y = unit(2, "lines"))
ggsave(plot = g, width = 12, height = 7, dpi = 300, filename = "tissue_example.pdf")		

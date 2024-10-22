setwd("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/Mdata/Example_figure")
library(ggpubr)		
library(RColorBrewer)	
library(ggplot2)
library(stringr)
theme_set(theme_classic())		
mpg<-read.table("total.txt",header=F,sep="\t")		
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
mpg$xx <- factor(mpg$xx, levels = c("AML","LIHC","COAD"))	#levels = rev(unique(mpg$xx)))		
mpg$zz <- factor(mpg$zz, levels = c('MCC','SCC','NCC'))	#rev(
mpg<-mpg[!is.na(mpg$zz),]	
mpg<-mpg[!is.na(mpg$yy),]	

g <- ggplot(mpg, aes(x=zz, y=yy, fill=zz)) + 
  # geom_boxplot(aes(middle = mean(yy))) +
  geom_boxplot(fatten = 2,outlier.shape = NA) +
  stat_summary(fun.y = mean, geom = "point", shape=20, size=3, color="black", fill="black",
               position = position_dodge(1.1))+
  facet_wrap(~xx, scale="free",nrow=1,labeller = as_labeller(c(
    "AML"="Acute myeloid leukemia","LIHC"="Liver hepatocellular carcinoma",
    "COAD"="Colon adenocarcinoma")))+
  scale_fill_manual(values = c("#1B9E77","#D95F02","#7570B3"))+
  xlab("") +	coord_cartesian(ylim = c(0, 0.4))+
  ylab("Methylation difference")+
  labs(fill="CpG type")+theme(legend.position = "none",
                              axis.title = element_text(size = 12),	strip.text = element_text(size = 14),
                              axis.text = element_text(size = 12),panel.spacing.x = unit(2, "lines"))
# +
#   stat_summary(fun.y="mean", geom="point", size=5,
#                position=position_dodge(width=0.75), color="white")
ggsave(plot = g, width = 12, height = 3.2, dpi = 300, filename = "disease_example.pdf")		

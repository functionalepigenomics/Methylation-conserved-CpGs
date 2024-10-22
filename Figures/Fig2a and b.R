setwd("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/Mdata/Monkey_probes_selection/Correlation_MC")
require('matrixStats')
datanew_5<-read.table("../../Old_Monkey_probes_selection/datanew_5.txt",header=T,row.names = 1)
test<-datanew_5
sd1<-c(NA,NA);sd2<-c(NA,NA);b<-c(NA,NA);c<-c(NA,NA);dif<-c(NA,NA)
for (i in 1:nrow(test)){
  x1<-mean(as.numeric(test[i,1:5]),na.rm=T)
  x2<-mean(as.numeric(test[i,6:11]),na.rm=T)
  x3<-mean(as.numeric(test[i,12:17]),na.rm=T)
  x4<-mean(as.numeric(test[i,18:23]),na.rm=T)
  x5<-mean(as.numeric(test[i,24:32]),na.rm=T)
  sd1[i]<-as.numeric(sd(c(x1,x2,x3,x4,x5)))
  b[i]<-max(mean(as.numeric(x1)),mean(as.numeric(x2)),mean(as.numeric(x3)),mean(as.numeric(x4)),mean(as.numeric(x5)))
  c[i]<-min(mean(as.numeric(x1)),mean(as.numeric(x2)),mean(as.numeric(x3)),mean(as.numeric(x4)),mean(as.numeric(x5)))
  dif[i]<-(b[i]-c[i])   ##### most important change
  sd2[i]<-as.numeric(sd(as.matrix(test[i,24:32]),na.rm = T))
}
dd<-data.frame(rownames(test),sd1,sd2)
rownames(dd)<-dd$rownames.test.

MC_CpG<-read.table("../../../MSCC.txt")
dd1<-dd[as.character(MC_CpG$V1),]
summary(dd1$sd1)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0005484 0.0158929 0.0243191 0.0250996 0.0336002 0.0499981 
summary(dd1$sd2)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001177 0.011185 0.016761 0.017355 0.021801 0.158882 

### Get SC-CpG
dd2<-dd[-which(rownames(dd) %in% as.character(MC_CpG$V1)),]
summary(dd2$sd1)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.004794 0.035361 0.074233 0.123221 0.162807 0.957145
summary(dd2$sd2)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0009371 0.0096873 0.0183220 0.0222573 0.0293540 0.3637236 

dd1$type<-"MCC"
dd2$type<-"SCC"

df<-rbind(dd1,dd2)
df$type <- factor(df$type, levels = c("MCC", "SCC"))

library(ggplot2)
require(ggExtra)
o4 <-  ggplot(df, aes(sd1, sd2)) + theme_classic()+
  geom_point(aes(colour = type),
             alpha = .01
             ) +labs(x="Interspecies\nmethylation variability (SD)",y="Interindividual\nmethylation variability (SD)")+
  scale_color_manual(values=c("#1B9E77","#D95F02"))+xlim(0, 0.4)+ylim(0, 0.3)+
  geom_smooth(method='lm', se=FALSE, color='grey')+
  theme(axis.ticks=element_line(colour = "black", size=1),
        panel.grid.major = element_blank(),axis.text=element_text(size=rel(1.2)),
        panel.grid.minor =element_blank(),axis.title=element_text(size=rel(1.2)), legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, size=1))
  # o5<-ggMarginal(
  # o4,
  # type = "boxplot",
  # margins = 'y',#'both',
  # size = 5,#varwidth = FALSE,notch = FALSE, notchwidth = 0.5,
  # groupColour = TRUE,
  # groupFill = TRUE) 
ggsave(plot = o4, width = 4, height = 3.5, dpi = 300, filename = "sd_correlation_type.pdf")

library(ggpubr)		
library(RColorBrewer)	
library(ggplot2)
library(stringr)
theme_set(theme_classic())	
g <- ggplot(df, aes(x=type, y=sd2, fill=type)) + 
  # geom_boxplot(aes(middle = mean(yy))) +
  geom_boxplot(fatten = NULL) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")+
  scale_fill_manual(values = c("#1B9E77","#D95F02"))+
  xlab("") +	
  ylab("Methylation variability (SD)")+
  labs(fill="CpG type")+theme(legend.position = "none",
                              axis.title = element_text(size = 12),	strip.text = element_text(size = 14),
                              axis.text = element_text(size = 12),panel.spacing.x = unit(2, "lines"))
ggsave(plot = g, width = 4, height = 3.5, dpi = 300, filename = "blood.pdf")
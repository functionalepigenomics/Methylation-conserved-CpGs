library(ggpubr)		
library(RColorBrewer)	
library(ggplot2)
library(stringr)
theme_set(theme_classic())	

df<-read.table("result_heredity.txt")
df$type<-NA
MSCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/MSCC.txt")
SCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/SCC.txt")
NCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/NCC.txt")
df$type[rownames(df) %in% as.character(MSCC$V1)]<-"MCC"
df$type[rownames(df) %in% as.character(SCC$V1)]<-"SCC"
df$type[rownames(df) %in% as.character(NCC$V1)]<-"NCC"
df<-df[!is.na(df$type),]
df$type <- factor(df$type, levels = c("MCC", "SCC","NCC"))
df$r2[df$r2<0]=0
df$r2[df$r2>1]=1

g <- ggplot(df, aes(x=type, y=r2, fill=type)) + 
  geom_boxplot(fatten = 2,outlier.shape = NA) +
  stat_summary(fun.y = mean, geom = "point", shape=20, size=3, color="black", 
               position = position_dodge(0.75))+
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))+
  xlab("") +	
  ylab("Heritability")+
  labs(fill="CpG type")+theme(legend.position = "none",
                              axis.title = element_text(size = 12),	strip.text = element_text(size = 14),
                              axis.text = element_text(size = 12),panel.spacing.x = unit(2, "lines"))
ggsave(plot = g, width = 4, height = 3.5, dpi = 300, filename = "heredity_boxplot.pdf")

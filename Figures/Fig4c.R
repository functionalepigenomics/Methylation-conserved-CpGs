library(ggplot2) 
library(plyr)
library(tidyverse)
df <-read.table("input",header=T,sep="\t")
df_sorted <- arrange(df, Cancer, type)
df_sorted$type <- factor(df_sorted$type, levels = c('MCC','SCC','NCC'))	#rev(
df_cumsum <- ddply(df_sorted, "Cancer",
                   transform, 
                   label_ypos=cumsum(value) - 0.5*value)
p <- ggplot(df_cumsum,aes(x=Cancer, y=value, fill=type)) +
  geom_bar(stat="identity", position='dodge')+
  geom_text(aes(label=value), vjust=1.6, position=position_dodge(width=0.9),
            color="white", size=3.5)+
  scale_fill_manual(values=c("#1B9E77","#7570B3", "#D95F02"))+
  theme_minimal()+ylim(0,50)+
  labs(y="Proportion of\ncancer-specific CpGs (%)")+
  theme(axis.line.y=element_line(size = 0.6),
        axis.line.x=element_blank(),
        axis.text.x=element_text(size=12,face="plain",margin = margin(r = 0.1)),
        axis.text.y=element_text(size=12,face="plain",margin = margin(r = 0.1)),
        axis.ticks.y=element_line(size = 1),
        #axis.ticks.length = unit(.001, "cm"),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12,face="plain",margin = margin(r = 0.1)),
        legend.position="right",
        legend.text = element_text(size=12,face="plain",margin = margin(r = 0.1)),
        legend.title = element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())#+
  #geom_hline(yintercept=16.72, color = "#ef8a62",size=1) # 16.72 #ef8a62   17.70 #67a9cf
ggsave(p,width = 7, height = 3, dpi = 300, filename = "bar_cancer.pdf")
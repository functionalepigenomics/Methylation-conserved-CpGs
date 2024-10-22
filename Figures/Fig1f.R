setwd("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/Mdata/Monkey_probes_selection")
require(ggpubr)
require(RColorBrewer)
dfm<-read.csv("input.gz",header=T,sep="\t")
dfm$cyl <- factor(dfm$cyl, levels = c("MCC", "SCC", "NCC"))
dfm$name <- factor(dfm$name, levels = c("Unmethylated", "Intermediately methylated", "Methylated"))
p<-ggplot(data=dfm, aes(x=cyl,y=mpg))  +theme_bw()+
  theme(#legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text = element_text(size=12),
        legend.title =element_blank()
  )+ylab("Proportion of CpGs")+
  #geom_rect(aes(),xmin =-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha = 0.1) +
  geom_bar(aes(fill = factor(name)),stat = "identity") +
  scale_fill_manual(values=c("#deebf7","#9ecae1","#3182bd"))#brewer.pal(3,"Dark2"))
# +
#   facet_wrap(~ name, nrow=1)
ggsave("bar.chart_DNAmpattern.pdf",p, width=5, height=2.2)#, units="in", scale=3)


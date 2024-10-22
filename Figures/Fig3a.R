setwd("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/Mdata/Genotype-Methylation/Figure/Proportion_meQTL/")
require(ggpubr)
dfm<-read.csv("input",header=T,sep="\t")
#dfm$mpg=dfm$mpg*100
# p<-ggdotchart(dfm, x = "name", y = "mpg",
#               #group = "cyl", 
#               color = "cyl",
#               dot.size = 4,
#               palette = "Set1",
#               rotate = TRUE,ylab="Ratio of Genetic-associated CpGs across diverse CpG types",
#               sorting = "descending",
#               ggtheme = theme_bw()
# )+
#   theme(legend.position = c(0.88,0.15),
#         axis.title.y=element_blank(),
#         legend.text = element_text(size=12),
#         axis.title.x = element_text(size=12),
#         axis.text = element_text(size=11),
#         legend.title =element_blank()
#   )+#geom_hline(yintercept=0, linetype="dashed", color = "black",show.legend=F)+
#   scale_y_continuous(breaks = sort(c(seq(0, 0.35,length.out=8), 1)))
# ggsave("lollipop.chart.pdf",p, width=2, height=1.3, units="in", scale=3)

dfm$cyl <- factor(dfm$cyl, levels = c("MCC", "SCC", "NCC"))
p<-ggplot(data=dfm, aes(x=cyl,y=mpg))  +theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=14), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text = element_text(size=14),
        legend.title =element_blank()
  )+ylab("Proportion of mQTL-CpGs")+
  #geom_rect(aes(),xmin =-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha = 0.1) +
  geom_bar(aes(fill = factor(cyl)),stat = "identity") +
  scale_fill_manual(values=brewer.pal(3,"Dark2"))+
  facet_wrap(~ name, nrow=2)
ggsave("bar.chart_NEW.pdf",p, width=3, height=2, units="in", scale=3)


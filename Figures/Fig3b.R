setwd("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/Mdata/Genotype-Methylation/Figure/Genetic_increase_meth")
require(ggpubr)
dfm<-read.csv("input.log",header=T,sep="\t")
dfm$cyl <- factor(dfm$cyl, levels = c("MCC", "SCC", "NCC"))
theme_set(theme_bw())
p<-ggplot(dfm, aes(x=name, y=mpg,fill=name)) + 
  geom_boxplot(outlier.shape = NA) +
  theme(legend.position = c(0.10,0.95),
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=14), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text = element_text(size=14),
        legend.title =element_blank()
  )+ylab("Methylation variability")+scale_y_continuous(limits = c(0, 0.3))+
  scale_fill_manual(values=c("#e41a1c","#377eb8"))+
  facet_wrap(~tissue)
ggsave("boxplot_GeneticinvreaseMeth_withoutCPGtype_NEW.pdf",p, width=3, height=2, units="in", scale=3)

p<-c(NA,NA);n=0
for (i in unique(dfm$tissue)){n=n+1;
a<-dfm$mpg[dfm$tissue==i & dfm$name=='mQTL-CpG']
b<-dfm$mpg[dfm$tissue==i & dfm$name == 'Other-CpG']
## 1000 permuatations
MC<-a
SC<-b
m=0
for (j in 1:1000){
  d<-SC[sample(length(SC), length(MC))];
  if (mean(d)>=mean(MC)){m=m+1}
}
p[n]=m/1000 
}

# > p
# [1] 0 0 0 0 0 0


a<-c(NA,NA);n=0;b<-c(NA,NA)
for (i in unique(dfm$tissue)){n=n+1;
a[i]<-mean(dfm$mpg[dfm$tissue==i & dfm$name=='mQTL-CpG'])
b[i]<-mean(dfm$mpg[dfm$tissue==i & dfm$name == 'Other-CpG'])
}
# Blood      Brain     Buccal       LCLs     Saliva       Skin 
# NA         NA 0.03604289 0.03629986 0.04603116 0.09003605 0.07125738 0.07790413 
# Blood      Brain     Buccal       LCLs     Saliva       Skin 
# NA         NA 0.01771082 0.02010377 0.02982284 0.06840894 0.04160865 0.05880629

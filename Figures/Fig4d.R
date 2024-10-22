library(ggprism)
coad<-read.table("COAD.txt",header=T)
esca<-read.table("ESCA.txt",header=T)
paad<-read.table("PAAD.txt",header=T)
sba<-read.table("SBA.txt",header=T)
coad$cancer<-"COAD"
esca$cancer<-"ESCA"
paad$cancer<-"PAAD"
sba$cancer<-"SBA"

df<-rbind(coad,esca,paad,sba)
MSCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/MSCC.txt")
SCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/SCC.txt")
NCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/NCC.txt")
df$con[df$rownames.test. %in% as.character(MSCC$V1)]<-1
df$con[df$rownames.test. %in% as.character(SCC$V1)]<-0
df$con[df$rownames.test. %in% as.character(NCC$V1)]<-2
df$diff1<-abs(df$diff1)

df$con[df$con==1]<-"MCC"
df$con[df$con==0]<-"SCC"
df$con[df$con==2]<-"NCC"
df<-df[which(df$diff1 > 0.03 & df$FDR < 0.05),]
df<-df[!is.na(df$con),]
df$con <- factor(df$con, levels = c('MCC','SCC','NCC'))	#rev(
require(ggplot2)

p<-ggplot(data = df, aes(x=cancer, y=diff1,fill=con)) + 
  geom_boxplot(fatten = 2,outlier.shape = NA) +
  stat_summary(fun.y = mean, geom = "point", shape=20, size=3, color="black", 
               position = position_dodge(0.75))+
  theme_classic()+
  theme(axis.title.x = element_blank(),axis.text.y = element_text(size = 12),
        axis.title.y=element_text(size=12,face="plain",margin = margin(r = 0.1)),
        legend.text = element_text(size=12,face="plain",margin = margin(r = 0.1)),
        axis.text.x = element_text(size = 12),
        legend.position = "right",
        axis.line.y=element_line(size = 0.6),
        axis.ticks.y=element_line(size = 1),
        #axis.ticks.length = unit(.001, "cm"),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank())+
  scale_fill_manual(values=c("#1B9E77", "#D95F02","#7570B3"))+
  labs(y="Methylation difference of\ncancer-specific CpGs")+scale_color_discrete(name="")+
  coord_cartesian(ylim = c(0, 0.6))
  #+rotate_x_text(angle = 45)#+
  # add_pvalue(df_p_val, 
  #            xmin = "xmin", 
  #            xmax = "xmax",
  #            label = "{p.adj.signif}",
  #            tip.length = 0)
ggsave(p,file="Cancer_variability_example.pdf",width = 7, height = 3, dpi = 300)


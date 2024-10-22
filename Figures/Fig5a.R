library(ggprism)
blca<-read.table("../../BLCA/Re.txt")
brca<-read.table("../../BRCA/Re.txt")
coad<-read.table("../../COAD/Re.txt")
esca<-read.table("../../ESCA/Re.txt")
hnsc<-read.table("../../HNSC/Re.txt")
kirc<-read.table("../../KIRC/Re.txt")
kirp<-read.table("../../KIRP/Re.txt")
lihc<-read.table("../../LIHC/Re.txt")
paad<-read.table("../../PAAD/Re.txt")
prad<-read.table("../../PRAD/Re.txt")
# read<-read.table("../../READ/Re.txt")
sarc<-read.table("../../SARC/Re.txt")
thca<-read.table("../../THCA/Re.txt")
ucec<-read.table("../../UCEC/Re.txt")



blca<-blca[blca$V5<0.05,]
blca$cancer<-"BLCA"
brca<-brca[brca$V5<0.05,]
brca$cancer<-"BRCA"
coad<-coad[coad$V5<0.05,]
coad$cancer<-"COAD"
esca<-esca[esca$V5<0.05,]
esca$cancer<-"ESCA"
hnsc<-hnsc[hnsc$V5<0.05,]
hnsc$cancer<-"HNSC"
kirc<-kirc[kirc$V5<0.05,]
kirc$cancer<-"KIRC"
kirp<-kirp[kirp$V5<0.05,]
kirp$cancer<-"KIRP"
lihc<-lihc[lihc$V5<0.05,]
lihc$cancer<-"LIHC"
paad<-paad[paad$V5<0.05,]
paad$cancer<-"PAAD"
prad<-prad[prad$V5<0.05,]
prad$cancer<-"PRAD"
# read<-read[read$V5<0.05,]
# read$cancer<-"READ"
sarc<-sarc[sarc$V5<0.05,]
sarc$cancer<-"SARC"
thca<-thca[thca$V5<0.05,]
thca$cancer<-"THCA"
ucec<-ucec[ucec$V5<0.05,]
ucec$cancer<-"UCEC"



df<-rbind(blca,brca,coad,esca,hnsc,kirc,kirp,lihc,paad,prad,sarc,thca,ucec)
df$con<-NA
names(df)[3]<-"diff1"
MSCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/MSCC.txt")
SCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/SCC.txt")
NCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/NCC.txt")
df$con[df$V1 %in% as.character(MSCC$V1)]<-1
df$con[df$V1 %in% as.character(SCC$V1)]<-0
df$con[df$V1 %in% as.character(NCC$V1)]<-2
df$diff1<-abs(as.numeric(df$diff1))

df$con[df$con==1]<-"MCC"
df$con[df$con==0]<-"SCC"
df$con[df$con==2]<-"NCC"
# df<-df[which(df$diff1 > 0.03 & df$FDR < 0.05),]
df<-df[!is.na(df$con),]
df$con <- factor(df$con, levels = c('MCC','SCC','NCC'))	#rev(
require(ggplot2)

p<-ggplot(data = df, aes(x=cancer, y=diff1)) + geom_boxplot(aes(fill=con),outlier.shape = NA)+theme_classic()+
  theme(axis.title.x = element_blank(),axis.text.y = element_text(size = 12),
        axis.title.y=element_text(size=12,face="plain",margin = margin(r = 0.1)),
        legend.text = element_text(size=12,face="plain",margin = margin(r = 0.1)),
        axis.text.x = element_text(size = 12),
        legend.position = "right",
        legend.title = element_blank(),
        axis.line.y=element_line(size = 0.6),
        axis.ticks.y=element_line(size = 1),
        #axis.ticks.length = unit(.001, "cm"),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank())+
  scale_fill_manual(values=c("#1B9E77", "#D95F02","#7570B3"))+
  labs(y="Regression coefficient")+scale_color_discrete(name="")+scale_y_continuous(limits = c(0, 1))
#+rotate_x_text(angle = 45)#+
# add_pvalue(df_p_val, 
#            xmin = "xmin", 
#            xmax = "xmax",
#            label = "{p.adj.signif}",
#            tip.length = 0)
ggsave(p,file="Cancer_variability_example.pdf",width = 12, height = 3, dpi = 300)


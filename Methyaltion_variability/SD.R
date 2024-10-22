setwd("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/Mdata/Human_population/Blood")
# library("gplots")
# #library("heatmap.plus")
# library("RColorBrewer")
# zz=gzfile('~/KoborLab/kobor_space/zdong/Monkey/Probes3/Mdata/Disorder/Cancer/CESC/datanew_5.txt.gz','rt')  
# test<-read.table(zz,header=T)
# rownames(test)<-test$TargetID
# test1=test[,grepl('Control',names(test))] 

# ##==========================Methylation Variance====================================
# require('matrixStats')
# test2<-test1
# ave1<-c(1,1);ave2<-c(1,1); ave3<-c(1,1);
# p1<-c(1,1);p2<-c(1,1);p3<-c(1,1);n=0;
# m<-rowSds(as.matrix(test2,na.rm=TRUE))
# m<-data.frame(rownames(test2),m)
# write.table(m,file="SD.txt",quote=F,col.names = F,row.names = F)
# # }


data<-read.table("SD.txt")
rownames(data)<-data$V1
colnames(data)[2]<-"diff1"
data$con<-NA
MSCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/MSCC.txt")
SCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/SCC.txt")
NCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/NCC.txt")
data$con[rownames(data) %in% as.character(MSCC$V1)]<-1
data$con[rownames(data) %in% as.character(SCC$V1)]<-0
data$con[rownames(data) %in% as.character(NCC$V1)]<-2
data<-data[!is.na(data$con),]
write.table(data,"../Tissue_figure/Blood",quote = F,sep="\t",col.names=F,row.names = F)
MC<-abs(data$diff1[data$con==1])
mean(MC,na.rm=TRUE) # 0.01869335
SC<-abs(data$diff1[data$con==0])
mean(SC,na.rm=TRUE) # 0.0234879
NC<-abs(data$diff1[data$con==2])
mean(NC,na.rm=TRUE) # 0.02478185

## 1000 permuatations
n=0
for (j in 1:1000){
  d<-SC[sample(length(SC), length(MC))];
  if (mean(d,na.rm=TRUE)<=mean(MC,na.rm=TRUE)){n=n+1}
}
p=n/1000 
p # 0

n1=0
for (j in 1:1000){
  d<-NC[sample(length(NC), length(MC))];
  if (mean(d,na.rm=TRUE)<=mean(MC,na.rm=TRUE)){n1=n1+1}
}
p1=n1/1000 
p1 # 0

n2=0
for (j in 1:1000){
  d<-NC[sample(length(NC), length(SC))];
  if (mean(d,na.rm=TRUE)<=mean(SC,na.rm=TRUE)){n2=n2+1}
}
p2=n2/1000 
p2 # 0

n2=0
for (j in 1:1000){
  d<-NC[sample(length(NC), length(SC))];
  if (mean(d,na.rm=TRUE)>=mean(SC,na.rm=TRUE)){n2=n2+1}
}
p2=n2/1000 
p2 # 1

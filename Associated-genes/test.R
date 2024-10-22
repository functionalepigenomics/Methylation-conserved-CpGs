require("matrixStats")
datanew_5<-read.table("~/KoborLab/kobor_space/zdong/Monkey/Probes3/Mdata/Disorder/Cancer/TCGA-BLCA/datanew_5.txt.gz",sep="\t",header = T)
datanew_5$con<-NA
test<-datanew_5
rownames(test)<-datanew_5$TargetID

# "TCGA.G4.6625.11A"

pheno<-read.table("README.txt",sep="\t",header = T)
pheno$bcr_sample_barcode<-gsub('-','\\.',  pheno$bcr_sample_barcode)
pheno<-pheno[pheno$bcr_sample_barcode %in% names(test),]


exp<-read.table(gzfile("TCGA.BLCA.sampleMap%2FHiSeqV2.gz"),header=T)
rownames(exp)<-exp$sample
names(exp)<-paste(names(exp), "A", sep="")
pheno<-pheno[pheno$bcr_sample_barcode %in% names(exp),]

test<-test[,names(test) %in% pheno$bcr_sample_barcode]
test[test>=1]<-0.9999
test[test<=0]<-0.0001
test=log2(test/(1-test))
exp<-exp[,names(exp) %in% pheno$bcr_sample_barcode]

pheno<-pheno[order(pheno$bcr_sample_barcode),]
test<-test[,order(names(test))]
exp<-exp[,order(names(exp))] # 394
pheno$type[grepl('.01[A-Z]$',pheno$bcr_sample_barcode)]<-1
pheno$type[grepl('.11[A-Z]$',pheno$bcr_sample_barcode)]<-0

pair<-read.table("cpg-gene.bed.gz",sep="\t")
pair<-pair[pair$V1 %in% rownames(test),]
pair<-pair[pair$V2 %in% rownames(exp),]

p<-c(NA,NA);beta<-c(NA,NA);
for (i in 1:nrow(pair)){
  df<-lm(as.numeric(test[rownames(test) %in% pair$V1[i],])~as.numeric(exp[rownames(exp) %in% pair$V2[i],])+
       pheno$age_at_initial_pathologic_diagnosis+pheno$gender+pheno$type)
  p[i]<-summary(df)$coefficients[2,4]
  beta[i]<-summary(df)$coefficients[2,1]
}
fdr<-p.adjust(p,method="fdr")
re<-data.frame(pair,beta,p,fdr)
write.table(re,file="Re.txt",sep="\t",row.names =F,col.names = F,quote=F)

#### ========== #####
re<-read.table("Re.txt")
names(re)[3]<-"beta"
names(re)[5]<-"fdr"
# re$fdr<-p.adjust(re$V4,method="bonferroni")
# write.table(re,file="Re.txt",sep="\t",row.names =F,col.names = F,quote=F)

MSCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/MSCC.txt")
SCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/SCC.txt")
NCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/NCC.txt")

m<-re[re$V1 %in% MSCC$V1,]
s<-re[re$V1 %in% SCC$V1,]
n<-re[re$V1 %in% NCC$V1,]

length(unique(m$V1)) # 7821
length(unique(s$V1)) # 92684
length(unique(n$V1)) # 149253

m1<-m[m$fdr<0.05,]
s1<-s[s$fdr<0.05,]
n1<-n[n$fdr<0.05,]

length(unique(m1$V1)) # 4487
length(unique(s1$V1)) # 53800
length(unique(n1$V1)) # 87415
length(unique(re$V1[re$fdr<0.05])) # 145702

length(unique(m1$V2)) # 3023
length(unique(s1$V2)) # 12619
length(unique(n1$V2)) # 15022
length(unique(re$V2[re$fdr<0.05])) # 15993
write.table(unique(m1$V2),file="mcc_genes.txt",sep="\t",row.names =F,col.names = F,quote=F)
write.table(unique(s1$V2),file="scc_genes.txt",sep="\t",row.names =F,col.names = F,quote=F)
write.table(unique(n1$V2),file="ncc_genes.txt",sep="\t",row.names =F,col.names = F,quote=F)

(length(unique(m1$V1))/length(unique(m$V1)))/(length(unique(s1$V1))/length(unique(s$V1))) # 0.9883625
(length(unique(m1$V1))/length(unique(m$V1)))/(length(unique(n1$V1))/length(unique(n$V1))) # 0.9795597
prop.test(x=c(length(unique(m1$V1)),length(unique(s1$V1))),n=c(length(unique(m$V1)),length(unique(s$V1)))) # 0.25
prop.test(x=c(length(unique(m1$V1)),length(unique(n1$V1))),n=c(length(unique(m$V1)),length(unique(n$V1)))) # 0.03727
prop.test(x=c(length(unique(s1$V1)),length(unique(n1$V1))),n=c(length(unique(s$V1)),length(unique(n$V1)))) # 0.01153
# > p.adjust(c(0.25,0.03727,0.01153),method="fdr")
# [1] 0.250000 0.055905 0.034590

mean(abs(m1$beta)) # 0.2345148
mean(abs(s1$beta)) # 0.2406799
mean(abs(n1$beta)) # 0.2463493

MC<-abs(m1$beta)
SC<-abs(s1$beta)
n=0
for (j in 1:1000){
  d<-SC[sample(length(SC), length(MC))];
  if (mean(d)<=mean(MC)){n=n+1}
}
p=n/1000 
p # 0.009

MC<-abs(m1$beta)
SC<-abs(n1$beta)
n=0
for (j in 1:1000){
  d<-SC[sample(length(SC), length(MC))];
  if (mean(d)<=mean(MC)){n=n+1}
}
p=n/1000 
p # 0

MC<-abs(s1$beta)
SC<-abs(n1$beta)
n=0
for (j in 1:1000){
  d<-SC[sample(length(SC), length(MC))];
  if (mean(d)<=mean(MC)){n=n+1}
}
p=n/1000 
p # 1

write.table(m1,file="mcc.txt",sep="\t",row.names =F,col.names = F,quote=F)
write.table(s1,file="scc.txt",sep="\t",row.names =F,col.names = F,quote=F)
write.table(n1,file="ncc.txt",sep="\t",row.names =F,col.names = F,quote=F)
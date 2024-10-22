setwd("~/KoborLab/kobor_space/kandy/home/zdong/Monkey/Probes3/Genotype-Methylation/Whole-blood/Ilmn.strand.flip/Correlation")

##================================Generate methylation input file======================================
a<-read.table("../../datanew_5.txt",header = T)
## before that you need to built README file by test.sh
# ethnicity.colname<-read.table("../../README.txt",header = T,sep = "\t")
# unmatch<-c(0,0)
# n=1;
# for (j in 1:length(colnames(a))){
#   for (i in 1:length(rownames(ethnicity.colname))){
#     if (colnames(a)[j] == ethnicity.colname$Sample[i]) {
#       colnames(a)[j]<-as.character(ethnicity.colname$Sample_ID[i])} #else {unmatch[n]<-colnames(test)[i+1];n=n+1}
#   }
# }

### get the list from SNP header
genotype<-read.table("Genotye_samplenames.txt",header=T,sep="\t")
colnames(genotype)[1]<-"TargetID"
## modified the column names in genotype file
# colnames(genotype)<-gsub(".GType", "",colnames(genotype))
# colnames(genotype)<-gsub(".R1", "",colnames(genotype))
b<-colnames(genotype)
c<-b[b%in%colnames(a)]
## remove duplicated samples
dup<-colnames(a[ , grepl( "_rep" , names(a) ) ])
c<-c[!c%in%dup]
## remove those samples in genotype data
Mismatch_snp<-b[!b%in%colnames(a)]
d<-a[,c]
write.table(d,file="Input_methylation.txt",quote = F,sep="\t",row.names = F)

##================================Generate genotype input file======================================
## merge into all.transpose.txt
#cat chr*.transpose.txt > all.transpose.txt
#sort -k2,2 -k3,3n all.transpose.txt > all.transpose_sorted.txt
#remove all colnames 
# cat Genotye_samplenames.txt all.transpose_sorted.txt > all.transpose.txt
genotype<-read.table("all.transpose.txt",header=T,sep="\t")
colnames(genotype)[1]<-"TargetID"
## modified the column names in genotype file
#colnames(genotype)<-gsub(".GType", "",colnames(genotype))
#colnames(genotype)<-gsub(".R1", "",colnames(genotype))
g1<-genotype[,c]
write.table(g1,file="Input_genotype_allchr.txt",quote = F,sep="\t",row.names = F)

##================================Generate covariates input file======================================
confounders<-read.table("../../a.log",header = T)
colnames(confounders)[1]<-"TargetID"
# unmatch<-c(0,0)
# n=1;
# for (j in 1:length(colnames(confounders))){
#   for (i in 1:length(rownames(ethnicity.colname))){
#     if (colnames(confounders)[j] == ethnicity.colname$Sample[i]) {
#       colnames(confounders)[j]<-as.character(ethnicity.colname$Sample_ID[i])} #else {unmatch[n]<-colnames(test)[i+1];n=n+1}
#   }
# }
confounders1<-confounders[,c]
confounder2<-confounders1[1:2,] ##!!!!! only keep gender and/or age into our analysis
write.table(confounder2,file="Input_covariate.txt",quote = F,sep="\t",row.names = F)





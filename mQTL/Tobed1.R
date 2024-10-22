a<-read.table("CpG_coordinate.txt")
## 10kb window for each side
a$V4=a$V3-1
a$V4[a$V4<1]<-1
a1<-a[,c("V1","V2","V4","V3")]
colnames(a1)<-c("geneid","chr","s1","s2")
write.table(a1,file="CpG_coordinate1.bed",col.names = T,row.names = F,quote = F,sep="\t")

b<-read.table("SNP_coordinate.txt")
b1<-b[,c("V3","V1","V2")]
b1$V1<-paste("chr", b1$V1, sep="")
colnames(b1)<-c("snp","chr","pos")
write.table(b1,file="SNP_coordinate1.bed",col.names =T,row.names = F,quote = F,sep="\t")
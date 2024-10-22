#setwd("~/KoborLab/kobor_space/zdong/Population_dataset/MHB_May3/Array_methylation")
library(wateRmelon)
library(dplyr) 
library(Rpdb)
library(sva)
data<-read.table("GSE55763/datanew_5.txt",header = T,row.names = 1,sep="\t")
data_dasen<-data[,2:(ncol(data)-1)]
# ## Check sex
# data1<-t(data)
# data1$estimatedSex <- sapply(1:nrow(data1), function(x){
#   XIC <- colMeans(data[c("cg03554089","cg12653510","cg05533223"),], na.rm = T)
#   if(XIC[x]>= 0.85){"M"} else{"F"}
# })
## Qunatile Normalization with Dasen
# data_dasen<-betaqn(data1)
# dim(data_dasen)
# # [1] 452728    133

test<-read.table("GSE73103/datanew_5.txt",header = T,row.names = 1,sep="\t")
test_dasen<-test[,2:(ncol(test)-1)]
# ## Check sex
# data1<-t(data)
# data1$estimatedSex <- sapply(1:nrow(data1), function(x){
#   XIC <- colMeans(data[c("cg03554089","cg12653510","cg05533223"),], na.rm = T)
#   if(XIC[x]>= 0.85){"M"} else{"F"}
# })
## Qunatile Normalization with Dasen
# test_dasen<-betaqn(test1)
# dim(test_dasen)
# # [1] 448892    192

train<-read.table("GSE77716/datanew_5.txt",header = T,row.names = 1, sep="\t")
rownames(train)<-train$TargetID
train_dasen<-train[,2:(ncol(train)-1)]
## Qunatile Normalization with Dasen
# train_dasen<-betaqn(train)
# dim(train_dasen)
# [1] 448892    192
#train_dasen<-train_dasen[,"GSM999376"]

valid<-read.table("David/datanew_5.txt",header = T,sep="\t")
rownames(valid)<-valid$TargetID
valid_dasen<-valid[,2:(ncol(valid)-1)]

common<-intersect(intersect(rownames(data_dasen),rownames(valid_dasen)),
                  intersect(rownames(test_dasen),rownames(train_dasen)))

data_dasen1<-data_dasen[as.character(common),]
test_dasen1<-test_dasen[as.character(common),]
train_dasen1<-train_dasen[as.character(common),]
valid_dasen1<-valid_dasen[as.character(common),]

total<-data.frame(data_dasen1,test_dasen1,train_dasen1,valid_dasen1)

subgroups<-c(rep("CEU",ncol(data_dasen1)))
batch<-c(rep("1",ncol(data_dasen1)))
pheno<-data.frame(subgroups,batch)
rownames(pheno)<-c(colnames(data_dasen1))

subgroups1<-c(rep("CEU",ncol(test_dasen1)))
batch1<-c(rep("2",ncol(test_dasen1)))
pheno1<-data.frame(subgroups1,batch1)
colnames(pheno1)<-c("subgroups","batch")
rownames(pheno1)<-c(colnames(test_dasen1))

subgroups2<-c(rep("CEU",ncol(train_dasen1)))
batch2<-c(rep("3",ncol(train_dasen1)))
pheno2<-data.frame(subgroups2,batch2)
colnames(pheno2)<-c("subgroups","batch")
rownames(pheno2)<-c(colnames(train_dasen1))

subgroups3<-c(rep("CEU",ncol(valid_dasen1)))
batch3<-c(rep("4",ncol(valid_dasen1)))
pheno3<-data.frame(subgroups3,batch3)
colnames(pheno3)<-c("subgroups","batch")
rownames(pheno3)<-c(colnames(valid_dasen1))

pheno4<-rbind(pheno,pheno1,pheno2,pheno3)
write.table(pheno4,file='pheno.txt',quote = F,sep='\t')
batch<-pheno4$batch
modcombat<-model.matrix(~1, data=pheno4)
combat_mydata= ComBat(dat=as.matrix(total), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
write.table(combat_mydata,file="combat_total.txt",quote=F,sep="\t")

require('matrixStats')
m<-rowSds(as.matrix(combat_mydata,na.rm=TRUE))
m<-data.frame(rownames(combat_mydata),m)
write.table(m,file="SD.txt",quote=F,col.names = F,row.names = F)





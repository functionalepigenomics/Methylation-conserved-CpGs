require(gganatogram)
library(gridExtra)
library(RColorBrewer)
male<-read.table("input_male.txt")
hgMale_key<-read.csv("male.txt",sep="\t")
test<-hgMale_key[hgMale_key$organ %in% male$V1[male$V2 == 1],]
male[which(male$V1 %in% hgMale_key$organ),] # 41
# write.table(male[which(male$V1 %in% hgMale_key$organ),],file="a.log",quote=F)
# test<-test[-which(test$organ %in% c("lung","spinal_cord","tonsil")),] ## change left ** to heart
gganatogram(data=test, fillOutline='white', organism='human', sex='male', fill="colour") +theme_void()

hgFemale_key<-read.csv("female.txt",sep="\t")
test1<-hgFemale_key[hgFemale_key$organ %in% male$V1[male$V2 == 0],]
# test1<-test1[-which(test1$organ %in% c("adipose_tissue","kidney","colon")),]
gganatogram(data=test1, fillOutline='white', organism='human', sex='female', fill="colour") +theme_void()
male[which(male$V1 %in% hgFemale_key$organ),] # 43
a<-rbind(test,test1)
length(unique(a$organ)) # 45
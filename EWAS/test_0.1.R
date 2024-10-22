setwd("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/Mcoad/Example_figure/Explain")
MSCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/MSCC.txt")
SCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/SCC.txt")
NCC<-read.table("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/NCC.txt")

#### ===== COAD =========
coad<-read.table("COAD.txt",header=T)
coad$con[coad$rownames.test. %in% as.character(MSCC$V1)]<-1
coad$con[coad$rownames.test. %in% as.character(SCC$V1)]<-0
coad$con[coad$rownames.test. %in% as.character(NCC$V1)]<-2
coad$diff1<-abs(coad$diff1)
nrow(coad[which(coad$diff1 > 0.10 & coad$FDR < 0.05 & coad$con == 1),])/nrow(coad[which(coad$con == 1),]) # 0.2469461
nrow(coad[which(coad$diff1 > 0.10 & coad$FDR < 0.05 & coad$con == 0),])/nrow(coad[which(coad$con == 0),]) # 0.2233543
nrow(coad[which(coad$diff1 > 0.10 & coad$FDR < 0.05 & coad$con == 2),])/nrow(coad[which(coad$con == 2),]) # 0.2297191
prop.test(x = c(nrow(coad[which(coad$diff1 > 0.10 & coad$FDR < 0.05 & coad$con == 1),]),nrow(coad[which(coad$diff1 > 0.10 & coad$FDR < 0.05 & coad$con == 0),])),
          n = c(nrow(coad[which(coad$con == 1),]),nrow(coad[which(coad$con == 0),])))
prop.test(x = c(nrow(coad[which(coad$diff1 > 0.10 & coad$FDR < 0.05 & coad$con == 1),]),nrow(coad[which(coad$diff1 > 0.10 & coad$FDR < 0.05 & coad$con == 2),])),
          n = c(nrow(coad[which(coad$con == 1),]),nrow(coad[which(coad$con == 2),])))
prop.test(x = c(nrow(coad[which(coad$diff1 > 0.10 & coad$FDR < 0.05 & coad$con == 0),]),nrow(coad[which(coad$diff1 > 0.10 & coad$FDR < 0.05 & coad$con == 2),])),
          n = c(nrow(coad[which(coad$con == 0),]),nrow(coad[which(coad$con == 2),])))

MC<-abs(coad$diff1[which(coad$diff1 > 0.10 & coad$FDR < 0.05 & coad$con == 1)])
mean(MC,na.rm=T) # 0.2142071
SC<-abs(coad$diff1[which(coad$diff1 > 0.10 & coad$FDR < 0.05 & coad$con == 0)])
mean(SC,na.rm=T) # 0.2092936
NC<-abs(coad$diff1[which(coad$diff1 > 0.10 & coad$FDR < 0.05 & coad$con == 2)])
mean(NC,na.rm=T) # 0.204144

## 1000 permuatations
n=0
for (j in 1:1000){
  d<-SC[sample(length(SC), length(MC))];
  if (mean(d,na.rm=T)<=mean(MC,na.rm=T)){n=n+1}
}
p=n/1000 
p # 0.999

n1=0
for (j in 1:1000){
  d<-NC[sample(length(NC), length(MC))];
  if (mean(d,na.rm=T)<=mean(MC,na.rm=T)){n1=n1+1}
}
p1=n1/1000 
p1 # 1

n2=0
for (j in 1:1000){
  d<-NC[sample(length(NC), length(SC))];
  if (mean(d,na.rm=T)<=mean(SC,na.rm=T)){n2=n2+1}
}
p2=n2/1000 
p2 # 1


### === ESCA =======
ESCA<-read.table("ESCA.txt",header=T)
ESCA$con[ESCA$rownames.test. %in% as.character(MSCC$V1)]<-1
ESCA$con[ESCA$rownames.test. %in% as.character(SCC$V1)]<-0
ESCA$con[ESCA$rownames.test. %in% as.character(NCC$V1)]<-2
ESCA$diff1<-abs(ESCA$diff1)
nrow(ESCA[which(ESCA$diff1 > 0.10 & ESCA$FDR < 0.05 & ESCA$con == 1),])/nrow(ESCA[which(ESCA$con == 1),]) # 0.1906417
nrow(ESCA[which(ESCA$diff1 > 0.10 & ESCA$FDR < 0.05 & ESCA$con == 0),])/nrow(ESCA[which(ESCA$con == 0),]) # 0.1850033
nrow(ESCA[which(ESCA$diff1 > 0.10 & ESCA$FDR < 0.05 & ESCA$con == 2),])/nrow(ESCA[which(ESCA$con == 2),]) # 0.1905763
prop.test(x = c(nrow(ESCA[which(ESCA$diff1 > 0.10 & ESCA$FDR < 0.05 & ESCA$con == 1),]),nrow(ESCA[which(ESCA$diff1 > 0.10 & ESCA$FDR < 0.05 & ESCA$con == 0),])),
          n = c(nrow(ESCA[which(ESCA$con == 1),]),nrow(ESCA[which(ESCA$con == 0),])))
prop.test(x = c(nrow(ESCA[which(ESCA$diff1 > 0.10 & ESCA$FDR < 0.05 & ESCA$con == 1),]),nrow(ESCA[which(ESCA$diff1 > 0.10 & ESCA$FDR < 0.05 & ESCA$con == 2),])),
          n = c(nrow(ESCA[which(ESCA$con == 1),]),nrow(ESCA[which(ESCA$con == 2),])))
prop.test(x = c(nrow(ESCA[which(ESCA$diff1 > 0.10 & ESCA$FDR < 0.05 & ESCA$con == 0),]),nrow(ESCA[which(ESCA$diff1 > 0.10 & ESCA$FDR < 0.05 & ESCA$con == 2),])),
          n = c(nrow(ESCA[which(ESCA$con == 0),]),nrow(ESCA[which(ESCA$con == 2),])))


MC<-abs(ESCA$diff1[which(ESCA$diff1 > 0.10 & ESCA$FDR < 0.05 & ESCA$con == 1)])
mean(MC,na.rm=T) # 0.1766307
SC<-abs(ESCA$diff1[which(ESCA$diff1 > 0.10 & ESCA$FDR < 0.05 & ESCA$con == 0)])
mean(SC,na.rm=T) # 0.180406
NC<-abs(ESCA$diff1[which(ESCA$diff1 > 0.10 & ESCA$FDR < 0.05 & ESCA$con == 2)])
mean(NC,na.rm=T) # 0.1736598

## 1000 permuatations
n=0
for (j in 1:1000){
  d<-SC[sample(length(SC), length(MC))];
  if (mean(d,na.rm=T)<=mean(MC,na.rm=T)){n=n+1}
}
p=n/1000 
p # 0.002

n1=0
for (j in 1:1000){
  d<-NC[sample(length(NC), length(MC))];
  if (mean(d,na.rm=T)<=mean(MC,na.rm=T)){n1=n1+1}
}
p1=n1/1000 
p1 # 0.985

n2=0
for (j in 1:1000){
  d<-NC[sample(length(NC), length(SC))];
  if (mean(d,na.rm=T)<=mean(SC,na.rm=T)){n2=n2+1}
}
p2=n2/1000 
p2 # 1


### === PAAD =====
PAAD<-read.table("PAAD.txt",header=T)
PAAD$con[PAAD$rownames.test. %in% as.character(MSCC$V1)]<-1
PAAD$con[PAAD$rownames.test. %in% as.character(SCC$V1)]<-0
PAAD$con[PAAD$rownames.test. %in% as.character(NCC$V1)]<-2
PAAD$diff1<-abs(PAAD$diff1)
nrow(PAAD[which(PAAD$diff1 > 0.10 & PAAD$FDR < 0.05 & PAAD$con == 1),])/nrow(PAAD[which(PAAD$con == 1),]) # 0.1372789
nrow(PAAD[which(PAAD$diff1 > 0.10 & PAAD$FDR < 0.05 & PAAD$con == 0),])/nrow(PAAD[which(PAAD$con == 0),]) # 0.1269851
nrow(PAAD[which(PAAD$diff1 > 0.10 & PAAD$FDR < 0.05 & PAAD$con == 2),])/nrow(PAAD[which(PAAD$con == 2),]) # 0.1138186
prop.test(x = c(nrow(PAAD[which(PAAD$diff1 > 0.10 & PAAD$FDR < 0.05 & PAAD$con == 1),]),nrow(PAAD[which(PAAD$diff1 > 0.10 & PAAD$FDR < 0.05 & PAAD$con == 0),])),
          n = c(nrow(PAAD[which(PAAD$con == 1),]),nrow(PAAD[which(PAAD$con == 0),])))
prop.test(x = c(nrow(PAAD[which(PAAD$diff1 > 0.10 & PAAD$FDR < 0.05 & PAAD$con == 1),]),nrow(PAAD[which(PAAD$diff1 > 0.10 & PAAD$FDR < 0.05 & PAAD$con == 2),])),
          n = c(nrow(PAAD[which(PAAD$con == 1),]),nrow(PAAD[which(PAAD$con == 2),])))
prop.test(x = c(nrow(PAAD[which(PAAD$diff1 > 0.10 & PAAD$FDR < 0.05 & PAAD$con == 0),]),nrow(PAAD[which(PAAD$diff1 > 0.10 & PAAD$FDR < 0.05 & PAAD$con == 2),])),
          n = c(nrow(PAAD[which(PAAD$con == 0),]),nrow(PAAD[which(PAAD$con == 2),])))

MC<-abs(PAAD$diff1[which(PAAD$diff1 > 0.10 & PAAD$FDR < 0.05 & PAAD$con == 1)])
mean(MC,na.rm=T) # 0.1646338
SC<-abs(PAAD$diff1[which(PAAD$diff1 > 0.10 & PAAD$FDR < 0.05 & PAAD$con == 0)])
mean(SC,na.rm=T) # 0.1686521
NC<-abs(PAAD$diff1[which(PAAD$diff1 > 0.10 & PAAD$FDR < 0.05 & PAAD$con == 2)])
mean(NC,na.rm=T) # 0.1672988

## 1000 permuatations
n=0
for (j in 1:1000){
  d<-SC[sample(length(SC), length(MC))];
  if (mean(d,na.rm=T)<=mean(MC,na.rm=T)){n=n+1}
}
p=n/1000 
p # 0.003

n1=0
for (j in 1:1000){
  d<-NC[sample(length(NC), length(MC))];
  if (mean(d,na.rm=T)<=mean(MC,na.rm=T)){n1=n1+1}
}
p1=n1/1000 
p1 # 0.017

n2=0
for (j in 1:1000){
  d<-NC[sample(length(NC), length(SC))];
  if (mean(d,na.rm=T)<=mean(SC,na.rm=T)){n2=n2+1}
}
p2=n2/1000 
p2 # 1

### === SBA =====
SBA<-read.table("SBA.txt",header=T)
SBA$con[SBA$rownames.test. %in% as.character(MSCC$V1)]<-1
SBA$con[SBA$rownames.test. %in% as.character(SCC$V1)]<-0
SBA$con[SBA$rownames.test. %in% as.character(NCC$V1)]<-2
SBA$diff1<-abs(SBA$diff1)
nrow(SBA[which(SBA$diff1 > 0.10 & SBA$FDR < 0.05 & SBA$con == 1),])/nrow(SBA[which(SBA$con == 1),]) # 0.1852335
nrow(SBA[which(SBA$diff1 > 0.10 & SBA$FDR < 0.05 & SBA$con == 0),])/nrow(SBA[which(SBA$con == 0),]) # 0.1637172
nrow(SBA[which(SBA$diff1 > 0.10 & SBA$FDR < 0.05 & SBA$con == 2),])/nrow(SBA[which(SBA$con == 2),]) # 0.1364418
prop.test(x = c(nrow(SBA[which(SBA$diff1 > 0.10 & SBA$FDR < 0.05 & SBA$con == 1),]),nrow(SBA[which(SBA$diff1 > 0.10 & SBA$FDR < 0.05 & SBA$con == 0),])),
          n = c(nrow(SBA[which(SBA$con == 1),]),nrow(SBA[which(SBA$con == 0),])))
prop.test(x = c(nrow(SBA[which(SBA$diff1 > 0.10 & SBA$FDR < 0.05 & SBA$con == 1),]),nrow(SBA[which(SBA$diff1 > 0.10 & SBA$FDR < 0.05 & SBA$con == 2),])),
          n = c(nrow(SBA[which(SBA$con == 1),]),nrow(SBA[which(SBA$con == 2),])))
prop.test(x = c(nrow(SBA[which(SBA$diff1 > 0.10 & SBA$FDR < 0.05 & SBA$con == 0),]),nrow(SBA[which(SBA$diff1 > 0.10 & SBA$FDR < 0.05 & SBA$con == 2),])),
          n = c(nrow(SBA[which(SBA$con == 0),]),nrow(SBA[which(SBA$con == 2),])))

MC<-abs(SBA$diff1[which(SBA$diff1 > 0.10 & SBA$FDR < 0.05 & SBA$con == 1)])
mean(MC,na.rm=T) # 0.2344033
SC<-abs(SBA$diff1[which(SBA$diff1 > 0.10 & SBA$FDR < 0.05 & SBA$con == 0)])
mean(SC,na.rm=T) # 0.2080661
NC<-abs(SBA$diff1[which(SBA$diff1 > 0.10 & SBA$FDR < 0.05 & SBA$con == 2)])
mean(NC,na.rm=T) # 0.1925915

## 1000 permuatations
n=0
for (j in 1:1000){
  d<-SC[sample(length(SC), length(MC))];
  if (mean(d,na.rm=T)<=mean(MC,na.rm=T)){n=n+1}
}
p=n/1000 
p # 1

n1=0
for (j in 1:1000){
  d<-NC[sample(length(NC), length(MC))];
  if (mean(d,na.rm=T)<=mean(MC,na.rm=T)){n1=n1+1}
}
p1=n1/1000 
p1 # 1

n2=0
for (j in 1:1000){
  d<-NC[sample(length(NC), length(SC))];
  if (mean(d,na.rm=T)<=mean(SC,na.rm=T)){n2=n2+1}
}
p2=n2/1000 
p2 # 1


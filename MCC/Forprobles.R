## Obtain the annotation data
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation.table = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

a<-annotation.table
a2<-a[a$Type=="II",]
a1<-a[a$Type=="I",]
a2$ProbeSeqB=a2$ProbeSeqA
a2$ProbeSeqA<-paste0(a2$ProbeSeqA,"A")
a2$ProbeSeqB<-paste0(a2$ProbeSeqB,"G")
Probes_Reverse<-rbind(a1,a2)
write.table(Probes_Reverse, file = "Probes.txt", append = FALSE, quote = F, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")


# a<-annotation.table[annotation.table$strand=="+",]
# a2<-a[a$Type=="II",]
# a1<-a[a$Type=="I",]
# a2$ProbeSeqB=a2$ProbeSeqA
# a2$ProbeSeqA<-paste0(a2$ProbeSeqA,"A")
# a2$ProbeSeqB<-paste0(a2$ProbeSeqB,"G")
# Probes_Forward<-rbind(a1,a2)
# write.table(Probes_Forward, file = "Probes_Forward.txt", append = FALSE, quote = F, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")

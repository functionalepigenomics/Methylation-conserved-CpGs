coad<-read.table("../../cancer-driver-gene/coad_cancergenes.txt")
esca<-read.table("../../cancer-driver-gene/esca_cancergenes.txt")
paad<-read.table("../../cancer-driver-gene/paad_cancergenes.txt")
# read<-read.table("../../cancer-driver-gene/read_cancergenes.txt")
lihc<-read.table("../../cancer-driver-gene/lihc_cancergenes.txt")
non<-read.table("../../cancer-driver-gene/non_cancergenes.txt")

overlap<-list(COAD=coad$V1,ESCA=esca$V1,
              PAAD=paad$V1,#READ=read$V1,
              LIHC=lihc$V1,NonGI=non$V1)
library("UpSetR")
library("ComplexHeatmap")
m = make_comb_mat(overlap)
pdf(file = "upset.pdf",width = 9, height = 3.4)
UpSet(m, top_annotation = upset_top_annotation(m, 
                                               annotation_name_rot = 90,
                                               annotation_name_side = "right",
                                               axis_param = list(side = "right")))
dev.off()
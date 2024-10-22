setwd("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/Mdata/Monkey_probes_selection")
p_sort<-read.table("../Old_Monkey_probes_selection/varibale.txt")
p_cutoff <- p_sort[1:length(p_sort$ID[p_sort$'q.meta' > 0.1]),]
## Filter the Methylation difference less than 0.1
p_sort2 <- p_cutoff[order(p_cutoff$Diff_beta),]
p_cutoff2 <- p_sort2[1:length(p_sort2$ID[p_sort2$Diff_beta < 0.05]),]
#p_cutoff2 <- p_cutoff2[order(p_cutoff2$`q.meta`,decreasing =T),]
#p_cutoff2<-p_cutoff
write.table(p_cutoff2,file="conservedcpgs.txt",quote=F)
rownames(p_sort)<-p_sort$ID
rownames(p_cutoff2)<-p_cutoff2$ID
SCC<-p_sort[-which(rownames(p_sort) %in% rownames(p_cutoff2)),]
write.table(rownames(p_cutoff2),file="../../MSCC.txt",quote=F,col.names = F,row.names = F)
write.table(rownames(SCC),file="../../SCC.txt",quote=F,col.names = F,row.names = F)
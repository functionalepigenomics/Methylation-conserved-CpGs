re<-read.table("Re.txt")
names(re)[3:5]<-c("dif","p","fdr")
re<-re[order(re$fdr,decreasing =T),]


### pca analysis
library(ggrepel)
# add a column of NAs
re$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
# re$diffexpressed[re$fdr < 0.05 & re$V1 %in% MSCC$V1] <- "MCC"
# # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
# re$diffexpressed[re$fdr < 0.05 & re$V1 %in% SCC$V1] <- "SCC"
# re$diffexpressed[re$fdr < 0.05 & re$V1 %in% NCC$V1] <- "NCC"
# re$diffexpressed[re$fdr < 0.05 & re$V1 %in% MSCC$V1] <- "MCC"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
re$diffexpressed[re$fdr < 0.05 & re$dif>0] <- "UP"
re$diffexpressed[re$fdr < 0.05 & re$dif<0] <- "DOWN"

re$diffexpressed <- factor(re$diffexpressed, levels = c('UP','DOWN','NO'))
g<-ggplot(data=re, aes(x=dif, y=-log10(p), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  #geom_text_repel() +
  scale_color_manual(values=c("#d7191c","#2b83ba", "#bdbdbd")) +
  # geom_vline(xintercept=c(-0.1, 0.1), col="black",linetype ="dashed") + 
  geom_hline(yintercept=-log10(head(re$p[re$fdr < 0.05],n=1)), col="black",linetype ="dashed")+
  theme(legend.position = "none",
        axis.text = element_text(size = 11)
  )+xlab("Regression coefficient")
ggsave(g,filename = "valcano.pdf",width = 4,height = 3.4)
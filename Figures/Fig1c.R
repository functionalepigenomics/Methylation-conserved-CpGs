setwd("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/Mdata/Monkey_probes_selection")
library("gplots")
#library("heatmap.plus")
library("RColorBrewer")
test<-read.table("../Old_Monkey_probes_selection/datanew_5.txt",header=T,row.names = 1)
con<-read.table("../../MSCC.txt")
test<-test[as.character(con$V1),]

coln<-read.table("../Old_Monkey_probes_selection/b.log",sep = "\t",header = T)
colnames(test)<-colnames(coln)
sum(is.na(test)) ## 441 NA values
## Imputation
row_average <- as.numeric(rowMeans(test[,1:5], na.rm=TRUE))
for (i in 1:length(row.names(test))) {
  test[i,1:5][is.na(test[i,1:5])] <- row_average[i]
}
row_average <- as.numeric(rowMeans(test[,6:11], na.rm=TRUE))
for (i in 1:length(row.names(test))) {
  test[i,6:11][is.na(test[i,6:11])] <- row_average[i]
}
row_average <- as.numeric(rowMeans(test[,12:17], na.rm=TRUE))
for (i in 1:length(row.names(test))) {
  test[i,12:17][is.na(test[i,12:17])] <- row_average[i]
}
row_average <- as.numeric(rowMeans(test[,18:23], na.rm=TRUE))
for (i in 1:length(row.names(test))) {
  test[i,18:23][is.na(test[i,18:23])] <- row_average[i]
}
row_average <- as.numeric(rowMeans(test[,24:32], na.rm=TRUE))
for (i in 1:length(row.names(test))) {
  test[i,24:32][is.na(test[i,24:32])] <- row_average[i]
}

condition_colors <- unlist(lapply(colnames(test),function(x){
  if(grepl('Chimpanzee',x)) '#E7298A'#'#FF9900' #pink
  else if(grepl('Bonobo',x)) '#66A61E'#'#800080' #grey
  else if(grepl('^G.',x)) '#E6AB02'#'#008000'
  else if(grepl('Human',x)) '#A6761D'#'#0000FF'
  else if(grepl('Orangutan',x)) '#666666'#'#FF0000'
}))
input<-as.matrix(test)
# myCol <- colorRampPalette(c("blue", "green", "yellow", "orange", "red"))(n = 499)
# myBreaks <- c(seq(0,0.2,length=100),  
#               seq(0.2,0.4,length=100),
#               seq(0.4,0.6,length=100),
#               seq(0.6,0.8,length=100),
#               seq(0.8,1,length=100))
my_palette <- colorRampPalette(c("yellow","orange", "red"))(n = 299)
col_breaks = c(seq(0,0.33,length=100),  # forestgreen
               seq(0.33,0.67,length=100), # yellow
               seq(0.67,1,length=100))    # red
pdf("test5.pdf")
heatmap.2(input,na.rm=T,Colv=NA,lhei = c(1,5),
          key.title=NA,key.ylab=NA,
          density="density",#tracecol="#303030",
          labRow=NA,labCol=NA,trace="none", #density="none", 
          col=my_palette,symbreaks = F,cexRow=1, cexCol=0.2, margins = c(1,1),
          key.xtickfun = function() {
            #breaks = pretty(parent.frame()$breaks)
            #breaks = breaks[c(1,2,3,length(breaks))]
            list(at = c(0,0.25,0.5,0.75,1)
                 #labels = breaks
            )},
          #Euclidean distance with Ward's linkage
          distfun=function(x) dist(x, method="euclidean"), hclustfun=function(x) hclust(x, method="ward.D2"),
          #1 minus Pearson correlation distance with average linkage
          #distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="average"),
          keysize=1.1,key.xlab="Methylation level",ColSideColors=condition_colors, scale="none")#,hclust=function(x) hclust(x,method="average"),distfun=function(x) as.dist((1-cor(x))/2))
legend(0.092,0.954,bty = "n",text.width=c(0,0.142,0.143,0.134,0.1495),legend=c("Chimpanzee","Bonobo","Gorilla","Orangutan","Human"),
       inset=c(-0.05,-0.1),text.col = "black",lty=0,xpd=T,horiz=T,cex=0.75,text.font=1)
dev.off()


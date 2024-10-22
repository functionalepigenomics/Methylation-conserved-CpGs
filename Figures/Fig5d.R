setwd("~/KoborLab/kobor_space/zdong/T2T/Sequence-based/Cancer/10X/Pathway")
## After get all the values then generate the plot
require(ggpubr)
require(stringr)
dfm<-read.csv("input.txt",header=T,sep='\t')


# a<-aggregate(dfm[, 2], list(dfm$name), mean)
# a<-a[order(a$x,decreasing = T),]
# b<-a[1:20,]
# dfm<-dfm[dfm$name %in% b$Group.1,]

# firstup <- function(x) {
#   substr(x, 1, 1) <- toupper(substr(x, 1, 1))
#   x
# }
# dfm$name<-firstup(dfm$name)
# dfm$group <- factor(dfm$group, levels = c("COAD","ESCA","PAAD","LIHC"))
# dfm$cyl <- factor(dfm$cyl, levels = c('MCC vs SCC','MCC vs NCC'))
p<-ggdotchart(dfm, x = "name", y = "mpg",
              color = "#E64B35",#shape ="group",
              dot.size = 4,
              #palette = "Set1",
              rotate = TRUE,ylab="Fold enrichment",
              sorting = "none",
              ggtheme = theme_bw()
)+
  # geom_point(aes(shape=factor(group)))+
  scale_shape(solid = TRUE)+
  # scale_shape_manual(values=c(19,1) , guide = "none")+
  theme(legend.position = "top",axis.title.x = element_text(size=12),
        axis.title.y=element_blank(),axis.text = element_text(size=12),
        legend.title =element_blank(),legend.text = element_text(size=12)
  )+
  # scale_y_continuous(breaks = sort(c(seq(0, 30,length.out=5), 1)))+
  # scale_x_discrete(labels=function(dfm) str_wrap(dfm,width = 30))+
  # scale_y_continuous(limits = c(1, 4))+
  scale_colour_manual(values=c("#E64B35"))+
  # scale_x_discrete(limits = rev(dfm$name))
ggsave("Enrichmen_pahtway_coloncancer.pdf",p, width=8, height=5.94)
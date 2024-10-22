setwd("~/KoborLab/kobor_space/zdong/Monkey/New_Probes/Mdata/Monkey_probes_selection/PhyloEpiGenomics")
library(PhyloEpiGenomics)
data<-read.table("../datanew_5.txt",sep="\t",row.names = 1,header = T)
Chimpanzee<-rowMeans(data[,1:5],na.rm=T)
Bonobo<-rowMeans(data[,6:11],na.rm=T)
Gorilla<-rowMeans(data[,12:17],na.rm=T)
Orangutan<-rowMeans(data[,18:23],na.rm=T)
Human<-rowMeans(data[,24:32],na.rm=T)
meth_fraction_aln<-data.frame(Chimpanzee,Bonobo,Gorilla,Orangutan,Human)
discretization=list(c(-0.01,0.2),c(0.2,0.4),c(0.4,0.6),c(0.6,0.8),c(0.8,1))
meth_states_aln=discretize(meth_fraction_aln,discretization)
my_noJump_model=make_evolutionary_model(meth_states_aln,model="noJump",nstates=length(discretization))



#use the same rate class "1" for all 5 edges of the unrooted tree
branches_to_evolutionary_rate_classes = rep(1, nrow(ml_best_nucl_tree$edge))
nrow(ml_best_nucl_tree$edge)
branches_to_evolutionary_rate_classes

#reconstructs a tree on methylation data expanding/contracting the form of
#nucleotide tree
ml_meth_tree_based_on_nucl = maximize_tree_log_likelihood_extended(
  meth_states_aln,
  tree = ml_best_nucl_tree,
  Q = my_noJump_model$Q,
  my_noJump_model$pi,
  branches_to_evolutionary_rate_classes = branches_to_evolutionary_rate_classes
)

#resulting meth tree (left bottom) has the same proportions as the nucl tree
#(right top)
kronoviz(
  list(ml_best_nucl_tree, ml_meth_tree_based_on_nucl$tree),
  type = "unrooted",
  lab4ut = "axial",
  rotate.tree = 270
)
add.scale.bar()
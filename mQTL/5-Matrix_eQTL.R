# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

setwd("~/KoborLab/kobor_space/kandy/home/zdong/Monkey/Probes3/Genotype-Methylation/Whole-blood/Ilmn.strand.flip/Correlation")
library(MatrixEQTL)

## Location of the package with the data files.
#base.dir = find.package('MatrixEQTL');
base.dir = '~/KoborLab/kobor_space/kandy/home/zdong/Monkey/Probes3/Genotype-Methylation/Whole-blood/Ilmn.strand.flip/Correlation';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(base.dir, "/Input_genotype_allchr.txt", sep="");
snps_location_file_name = paste(base.dir, "/SNP_coordinate.sorted1.bed", sep="");

# Gene expression file name
expression_file_name = paste(base.dir, "/Input_methylation.txt", sep="");
gene_location_file_name = paste(base.dir, "/CpG_coordinate.sorted1.bed", sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste(base.dir, "/Input_covariate.txt", sep="");

# Output file name
output_file_name_cis = "result_cis.txt";
output_file_name_tra = "result_tra.txt";

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-2;
pvOutputThreshold_tra = 5e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e4;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = 0,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

# unlink(output_file_name_tra);
# unlink(output_file_name_cis);

## Results:

# cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
# cat('Detected local eQTLs:', '\n');
# show(me$cis$eqtls)
# cat('Detected distant eQTLs:', '\n');
# show(me$trans$eqtls)

re<-read.table("result_cis.txt",header = T)
re1<-re[which(re$FDR<0.05),]
#re2<-re1[unique(re1$gene),]
#write.table(re2,file="result_cis_uniqueCpG.txt",quote = F,sep="\t",row.names = F)
MC_CpG<-read.table("~/KoborLab/kobor_space/kandy/home/zdong/Monkey/Probes3/Mdata/Monkey_probes_selection/conservedcpgsprobes.txt")
SC_CpG<-read.table("~/KoborLab/kobor_space/kandy/home/zdong/Monkey/Probes3/Mdata/Monkey_probes_selection/nonconservedcpgprobes.txt")
NC_CpG<-read.table("~/KoborLab/kobor_space/kandy/home/zdong/Monkey/Probes3/Mdata/otherprobes.txt")
n_MC_CpG<-length(unique(re1$gene[re1$gene%in%MC_CpG$V1])) #973
n_MC_CpG
length(re1$gene[re1$gene%in%MC_CpG$V1])/length(unique(re1$gene[re1$gene%in%MC_CpG$V1])) #12.89517
write.table(re1$gene[re1$gene%in%MC_CpG$V1],file="result_cis_in_MC-CpG.txt",quote = F,sep="\t",row.names = F,col.names = F)
write.table(unique(re1$gene[re1$gene%in%MC_CpG$V1]),file="result_cis_uniqueCpG_in_MC-CpG.txt",quote = F,sep="\t",row.names = F,col.names = F)
n_SC_CpG<-length(unique(re1$gene[re1$gene%in%SC_CpG$V1])) #18078
n_SC_CpG
length(re1$gene[re1$gene%in%SC_CpG$V1])/length(unique(re1$gene[re1$gene%in%SC_CpG$V1])) #17.80828
write.table(re1$gene[re1$gene%in%SC_CpG$V1],file="result_cis_in_SC-CpG.txt",quote = F,sep="\t",row.names = F,col.names = F)
write.table(unique(re1$gene[re1$gene%in%SC_CpG$V1]),file="result_cis_uniqueCpG_in_SC-CpG.txt",quote = F,sep="\t",row.names = F,col.names = F)
n_NC_CpG<-length(unique(re1$gene[re1$gene%in%NC_CpG$V1])) #133214
n_NC_CpG
length(re1$gene[re1$gene%in%NC_CpG$V1])/length(unique(re1$gene[re1$gene%in%NC_CpG$V1])) #18.43136
write.table(re1$gene[re1$gene%in%NC_CpG$V1],file="result_cis_in_NC-CpG.txt",quote = F,sep="\t",row.names = F,col.names = F)
write.table(unique(re1$gene[re1$gene%in%NC_CpG$V1]),file="result_cis_uniqueCpG_in_NC-CpG.txt",quote = F,sep="\t",row.names = F,col.names = F)
QC_filter_CpG<-read.table("CpG_coordinate.sorted1.bed",header=T,row.names = 1)
ratio_MC<-n_MC_CpG/length(MC_CpG$V1[MC_CpG$V1%in%rownames(QC_filter_CpG)]) # 0.3094784
ratio_MC
ratio_SC<-n_SC_CpG/length(SC_CpG$V1[SC_CpG$V1%in%rownames(QC_filter_CpG)]) # 0.3829355
ratio_SC
ratio_NC<-n_NC_CpG/length(NC_CpG$V1[NC_CpG$V1%in%rownames(QC_filter_CpG)]) # 0.385542
ratio_NC
write.table(unique(re1$SNP[re1$gene%in%MC_CpG$V1]),file="result_cis_uniquemQTL_in_MC-mQTL.txt",quote = F,sep="\t",row.names = F,col.names = F)
write.table(unique(re1$SNP[re1$gene%in%SC_CpG$V1]),file="result_cis_uniquemQTL_in_SC-mQTL.txt",quote = F,sep="\t",row.names = F,col.names = F)
write.table(unique(re1$SNP[re1$gene%in%NC_CpG$V1]),file="result_cis_uniquemQTL_in_NC-mQTL.txt",quote = F,sep="\t",row.names = F,col.names = F)


###=====================================Tag-SNP based analysis====================================
re<-read.table("result_cis.txt",header = T)
re1<-re[which(re$FDR<0.05),]
tagsnp<-read.table('../Imputation_Minimac3andeagle_hrc.r1.1.2016/Tag-snp_selection/all.snp.purifed.imputed.nonduplicate.txt')
re1<-re1[re1$SNP %in% tagsnp$V1,]
#re1<-re1[unique(re1$gene),]
#write.table(re1,file="result_cis_uniqueCpG.txt",quote = F,sep="\t",row.names = F)
MC_CpG<-read.table("~/KoborLab/kobor_space/kandy/home/zdong/Monkey/Probes3/Mdata/Monkey_probes_selection/conservedcpgsprobes.txt")
SC_CpG<-read.table("~/KoborLab/kobor_space/kandy/home/zdong/Monkey/Probes3/Mdata/Monkey_probes_selection/nonconservedcpgprobes.txt")
NC_CpG<-read.table("~/KoborLab/kobor_space/kandy/home/zdong/Monkey/Probes3/Mdata/otherprobes.txt")
length(re1$gene[re1$gene%in%MC_CpG$V1])/length(unique(re1$gene[re1$gene%in%MC_CpG$V1])) #3.218335
length(re1$gene[re1$gene%in%SC_CpG$V1])/length(unique(re1$gene[re1$gene%in%SC_CpG$V1])) #4.61488
length(re1$gene[re1$gene%in%NC_CpG$V1])/length(unique(re1$gene[re1$gene%in%NC_CpG$V1])) #4.660308
write.table(re1$gene[re1$gene%in%MC_CpG$V1],file="result_cis_tag-snpselection_in_MC-CpG.txt",quote = F,sep="\t",row.names = F,col.names = F)
write.table(re1$gene[re1$gene%in%SC_CpG$V1],file="result_cis_tag-snpselection_in_SC-CpG.txt",quote = F,sep="\t",row.names = F,col.names = F)
write.table(re1$gene[re1$gene%in%NC_CpG$V1],file="result_cis_tag-snpselection_in_NC-CpG.txt",quote = F,sep="\t",row.names = F,col.names = F)
system("bash count.sh",intern = TRUE)









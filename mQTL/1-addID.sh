for chr in $(seq 1 22); do
cut -f7- -d" " ../Imputation_Minimac3andeagle_hrc.r1.1.2016/Imputed.convert.chr$chr.raw > a.txt

sed -i "s/ /\t/g" a.txt
### need to check the number of colname in colname.txt and in Imputed.convert.chr$chr.raw files
#head -1 ../GSE24274_series_matrix.txt > colname.txt

## In this dataset, the colname file has been modified to match the nimber

awk -f ~/KoborLab/kobor_space/kandy/home/zdong/Monkey/Probes3/Genotype-Methylation/Brain_GSE112525/A1-b37.Ilmn_strand_flip/Correlation/transpose1.awk a.txt > chr$chr.transpose.txt

#cat colname.txt b.txt > chr$chr.transpose.txt

rm -f a.txt
done

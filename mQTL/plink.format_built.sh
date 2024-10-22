#!/bin/sh
## match genotype data and corresonding allale data
filename="../SNP.txt"
filenamestrand="../GSA-24v1-0_A2-b37.Ilmn.strand"
filename2="../GSA-24v1-0_A2.update_alleles.txt"
#sed 1d $filename | sed 's/-.*[0-9]//' | sort -k1,1 > c.log
sort -k1,1 $filename > c.log
sort -k1,1 $filename2 > d.log
filename="c.log"
filename2="d.log"
awk '{print $1}' c.log > a.log 
awk '{print $1}' d.log > b.log

## After match
if diff a.log b.log >& /dev/null;
 then
 # MAP FORMAT:
 # chromosome (1-22, X, Y or 0 if unplaced)
 # rs# or snp identifier
 # Genetic distance (morgans)
 # Base-pair position (bp units)
 # Generate MAP file
 awk '{print ""$1"\t"$2"\t"$1"\t0\t"$3""}' $filenamestrand  > plink.map

 ## pay attention to SNP132 .... 
 # make sure that .map are match with genotype data
 cut -f1 d.log > Sample.log
 awk 'FNR==NR{a[$1];next}($1 in a){print}' Sample.log plink.map  > plink1.map

 cut -f1 plink1.map > Sample.log
 awk 'FNR==NR{a[$1];next}($1 in a){print}' Sample.log c.log > c_matched.log
 awk 'FNR==NR{a[$1];next}($1 in a){print}' Sample.log d.log > d_matched.log
 sort -k1,1 plink1.map  > plink_sorted.map 
 awk '{print $1}' plink_sorted.map > Sample_a.log
 awk '{print $1}' d_matched.log > Sample_b.log
 # Generate the matched MAP file and PED file
 if diff Sample_a.log Sample_b.log >& /dev/null;
  then cut -f2- plink_sorted.map > plink_matched.map
  cut -f1,3 d_matched.log | sed 's/ /\t/g' > d_new.log  
  paste d_new.log c_matched.log | cut -f1-3,5- | sed 's/AA/0\t/g' | sed 's/AB/1\t/g' | sed 's/BB/2\t/g' | sed 's/NA/3\t/g'> myRaw.txt  ## only for this data, missing genotype was shown as NA, others as NC
  Rscript ped_built.R
 fi
fi

# To convert .ped and .map in Plink binary format and pre-quality check
plink --file plink_matched --maf 0.05 \
	--mind 0.1 \
	--geno 0.1 \
	--hwe 0.001 \
#	--exclude mysnps.txt \
	--make-bed --out plink_matched_binary

# split by chromosome
for chr in $(seq 1 22); do
     plink --bfile plink_matched_binary \
           --chr $chr \
           --make-bed --out plink_matched_binary_chr$chr;
     ~/KoborLab/kobor_space/kandy/home/zdong/Monkey/Probes3/Genotype-Methylation/LCLs_GSE24277/updata_build.sh plink_matched_binary_chr$chr $filenamestrand plink_matched_binary_updata_build_chr$chr
     rm plink_matched_binary_chr$chr.bed
     rm plink_matched_binary_chr$chr.bim
     rm plink_matched_binary_chr$chr.fam
done


# #####=============================================Method 1============================================
# #####################################
# ## check data ##
# #####################################
# ## This part come from https://bitbucket.org/sics-sb/gusto-imputation/src/e40b9d16219cb5b2659071da84c4538fcb6d78fa/4Phasing/familyPhasing.sh?at=master&fileviewer=file-view-default
# ## SHAPEIT automaltically checks the input genotype data prior to phasing. It throws warnings when it detects:
# ## Individuals or sites with a rate of missing data higher than 5%.
# ## Monomorphic or singleton SNPs since they are not informative for phasing
# ## And errors when an individual or a site has a missing data rate of 100%. 
# for chr in $(seq 1 22); do
#     shapeit -check -B plink_matched_binary_updata_build_chr$chr --output-log plink_matched_binary_updata_build_chr$chr.checks
# done
# ## Note if bad probes and individuals were included, please remove them in phasing progress
# 
# ## Phase (with duohmm for autosomes) 
# ## 1000 genome v3 was used in prephasing
# for chr in $(seq 1 22); do
#     shapeit -check  \
#         -B plink_matched_binary_updata_build_chr$chr \
# 	-M ~/KoborLab/kobor_space/kandy/home/zdong/Sofwares/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/Geneticmap_1000v3/1000GP_Phase3/genetic_map_chr"$chr"_combined_b37.txt \
# 	##--exclude-snp gwas.subset.site \   ## The data is from plink_matched_binary_updata_build_chr$chr.checks.snp.mm and only positions were used 
# 	##--exclude-ind gwas.subset.inds \   ## The data is from plink_matched_binary_updata_build_chr$chr.checks.ind.mm and only individual IDs were used
# 	--input-ref ~/KoborLab/kobor_space/kandy/home/zdong/Sofwares/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/Reference_panel/1000GP_Phase3/1000GP_Phase3_chr$chr.hap.gz \
# 	~/KoborLab/kobor_space/kandy/home/zdong/Sofwares/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/Reference_panel/1000GP_Phase3/1000GP_Phase3_chr$chr.legend.gz \
# 	~/KoborLab/kobor_space/kandy/home/zdong/Sofwares/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/Reference_panel/1000GP_Phase3/1000GP_Phase3.sample \
# 	-W 5 \
# 	-O plink_matched_binary_prephased_chr$chr --states 200 --duohmm
# done
# 
# ## Convert to VCF
# for chr in $(seq 1 22); do
#     shapeit -convert \
#         --input-haps plink_matched_binary_prephased_chr$chr \
#         --output-vcf plink_matched_binary_prephased_chr$chr.vcf
#     vcf-sort plink_matched_binary_prephased_chr$chr.vcf | bgzip -c > plink_matched_binary_prephased_chr$chr.vcf.gz
# done
# ## Note: After phased, those vcfs can be imputed with Michigan Imputation Server (see below) directly using Minimac3 with Reference Panel (1000 genome v3 or HRC r1.1)


#####=============================================Method 2============================================
# prepare the input format for Michigan Imputation Server (https://imputationserver.sph.umich.edu/start.html#!pages/home) and where Eagle for phasing; Minimac3 for imputation using Reference Panel "HRC r1.1"
# Das S, Forer L, Schönherr S, et al. Next-generation genotype imputation service and methods. Nat Genet. 2016;48(10):1284-1287.
for chr in $(seq 1 22); do
     plink --bfile plink_matched_binary_updata_build_chr$chr --snps-only just-acgt --flip snp.strand.flip.txt --recode vcf --out plink_matched_binary_updata_build_chr$chr
     vcf-sort plink_matched_binary_updata_build_chr$chr.vcf | bgzip -c > plink_matched_binary_updata_build_chr$chr.vcf.gz
#     ~/KoborLab/kobor_space/kandy/home/zdong/checkVCF/checkVCF/checkVCF.py -r ~/KoborLab/kobor_space/kandy/home/zdong/checkVCF/checkVCF/hs37d5.fa -o out plink_matched_binary_updata_build_chr$chr.vcf.gz
     rm plink_matched_binary_updata_build_chr$chr.fam
     rm plink_matched_binary_updata_build_chr$chr.bed
     rm plink_matched_binary_updata_build_chr$chr.bim
     rm plink_matched_binary_updata_build_chr$chr.vcf
done
# After check the srand, do pahse and imputation now!!!!

## If there still include some strand-flip sites, just remove them in QC process:
     #--exclude mysnps.txt

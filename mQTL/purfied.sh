## obtain the target sample IDs
cut -f1,2 -d " " ../Imputed.chr22.ped > a.log
tr -s '\t'  '\n' < ../../Correlation/Genotye_samplenames.txt | tail -n +2 > b.log
paste b.log a.log > c.log
head -1 ../../Correlation/Input_genotype_allchr.txt | tr -s '\t'  '\n'  | tail -n +2 > d.log
## need modify c.log sample IDs in Brain data
awk 'FNR==NR{a[$1];next}($1 in a){print}' d.log c.log > e.log
cut -f2- e.log > keep_indviduals.txt
rm *.log

## 
for chr in $(seq 1 22); do
plink --file ../Imputed.chr$chr  --keep keep_indviduals.txt --make-bed --out Imputed.chr$chr
done
# ls Imputed.chr*.bed > bed.txt
# ls Imputed.chr*.bim > bim.txt
# ls Imputed.chr*.fam > fam.txt
# 
# paste bed.txt bim.txt fam.txt | column -s $'\t' -t > merge.temp.list
# 
# grep -v -w "chr1." merge.temp.list > merge.list
# remove Imputed.chr1. file in merge.temp.list and generat merge.list
# plink --bfile Imputed.chr1 --merge-list merge.list --make-bed --out all

for chr in $(seq 1 22); do
plink --bfile Imputed.chr$chr --indep-pairwise 50 5 0.8 --out tmp1
plink --bfile Imputed.chr$chr --extract tmp1.prune.in --recode --out purifed.imputed.dup.chr$chr
plink --bfile Imputed.chr$chr --extract tmp1.prune.in --recodeA --out snp.purifed.imputed.dup.chr$chr
head -1 snp.purifed.imputed.dup.chr$chr.raw |tr -s ' '  '\n'| tail -n +7 > snp.purifed.imputed.dup.chr$chr.log
### !!! remove duplicated ID
#plink --file purifed.imputed.dup.chr$chr --list-duplicate-vars ids-only suppress-first 
#plink --file purifed.imputed.dup.chr$chr -exclude plink.dupvar --recode --out purifed.imputed.chr$chr
### !!! remove duplicated genomic coordinate
#uniq -d purifed.imputed.dup.chr$chr.map > tmp2
#plink --file purifed.imputed.dup.chr$chr --exclude tmp2 --recode --out purifed.imputed.chr$chr
done
cat purifed.imputed.dup.chr*.map > all.purfied.imputated.snp.txt
sort all.purfied.imputated.snp.txt | uniq > all.purfied.imputated.snp.nonduplicate.txt
cat snp.purifed.imputed.dup.chr*.log > all.snp.purifed.imputed.dup.txt
sort all.snp.purifed.imputed.dup.txt | uniq > all.snp.purifed.imputed.nonduplicate.txt

rm *chr*

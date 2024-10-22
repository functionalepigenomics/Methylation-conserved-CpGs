for chr in $(seq 1 22); do
## unzip files with password
 unzip -P "R8-BLE8ytXjop3" chr_$chr.zip 
 gunzip chr$chr.dose.vcf.gz
## remove MAF < 0.05, exclude individuals with more than 10% missing genotypes,
## include only SNPs with a 90% genotyping rate (10% missing)
 plink --vcf chr$chr.dose.vcf --maf 0.05 --mind 0.1 --geno 0.1 --hwe 0.001 --recode --out Imputed.chr$chr
done

## convert ped format into the 012 format (o hom ancesteral), 1 het, 2 dom derived
for chr in $(seq 1 22); do
plink --file Imputed.chr$chr  --recodeA --out Imputed.convert.chr$chr
done

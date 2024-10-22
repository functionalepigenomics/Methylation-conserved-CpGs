#ls ../*.ped > ped.txt
#ls ../*.map > map.txt
#paste ped.txt map.txt | column -s $'\t' -t > merge.temp.list
#grep -v -w "chr1//." merge.temp.list > merge.list
for chr in $(seq 1 22); do
plink --file ../Imputed.chr$chr  --hardy  --out $chr
#remove header, substitute "/" to "tabs", calculate frequencies, output to new file
sed 1d $chr.hwe | \
sed 's:/:\t:g' | \
awk '{OFS="\t";print $1,$2,$3,$4,$5,(($6+$7/2)/($6+$7+$8)),(($8+$7/2)/($6+$7+$8)),$9,$10,$11}' \
> $chr.hwe.freq
done
cat *.hwe.freq > Hwe.freq.txt

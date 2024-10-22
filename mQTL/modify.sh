grep "Strand switch" statistics.txt > a.txt
#grep "Allele switch:" statistics.txt > b.txt
sed -i "s/FILTER - Strand switch and Allele switch: //g" a.txt
#sed -i "s/INFO - Allele switch: //g" b.txt
sed -i "s/FILTER - Strand switch: //g" a.txt
#cut -f1 b.txt > mysnps.txt
cut -f1 a.txt > snp.strand.flip.txt
#rm a.txt b.txt

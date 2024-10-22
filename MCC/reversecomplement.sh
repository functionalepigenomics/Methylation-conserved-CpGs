#perl -pe '/^>/ ? print "\n" : chomp' GCF_000364345.1_Macaca_fascicularis_5.0_genomic.fna | while read -r line; do if [[ "$line" == ">"* ]]; then echo $line; else tr "[ATGCatgcNn]" "[TACGtacgNn]"; fi; done | awk '{ if ($0~/^>/) { printf("%s",$0); } else { for(i=length;i!=0;i--) x=(x substr($0,i,1)) } { printf("%s\n",x); x=""; } }' > reverse_monkey.fasta 
sed -i -e 's/A/W/g' genome_mfa.GA_conversion_reversecomplement.fa
sed -i -e 's/T/A/g' genome_mfa.GA_conversion_reversecomplement.fa
sed -i -e 's/W/T/g' genome_mfa.GA_conversion_reversecomplement.fa
sed -i -e 's/C/W/g' genome_mfa.GA_conversion_reversecomplement.fa
sed -i -e 's/G/C/g' genome_mfa.GA_conversion_reversecomplement.fa
sed -i -e 's/W/G/g' genome_mfa.GA_conversion_reversecomplement.fa
sed -i -e 's/NG_/NC_/g' genome_mfa.GA_conversion_reversecomplement.fa
sed -i -e 's/CRGh37/GRCh37/g' genome_mfa.GA_conversion_reversecomplement.fa
sed -i -e 's/MA/MT/g' genome_mfa.GA_conversion_reversecomplement.fa
sed -i -e 's/CL/GL/g' genome_mfa.GA_conversion_reversecomplement.fa
sed -i -e 's/Gomplete/Complete/g' genome_mfa.GA_conversion_reversecomplement.fa


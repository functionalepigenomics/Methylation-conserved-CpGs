#Fields are space-separated, and all 21 are required for PSL File Format

# 1.matches - Number of matching bases that aren't repeats.
# 2.misMatches - Number of bases that don't match.
# 3.repMatches - Number of matching bases that are part of repeats.
# 4.nCount - Number of 'N' bases.
# 5.qNumInsert - Number of inserts in query.
# 6.qBaseInsert - Number of bases inserted into query.
# 7.tNumInsert - Number of inserts in target.
# 8.tBaseInsert - Number of bases inserted into target.
# 9.strand - defined as + (forward) or - (reverse) for query strand. In mouse, a second '+' or '-' indecates genomic strand.
# 10.qName - Query sequence name.
# 11.qSize - Query sequence size.
# 12.qStart - Alignment start position in query.
# 13.qEnd - Alignment end position in query.
# 14.tName - Target sequence name.
# 15.tSize - Target sequence size.
# 16.tStart - Alignment start position in target.
# 17.tEnd - Alignment end position in target.
# 18.blockCount - Number of blocks in the alignment.
# 19.blockSizes - Comma-separated list of sizes of each block.
# 20.qStarts - Comma-separated list of start position of each block in query.
# 21.tStarts - Comma-separated list of start position of each block in target.

#Filter.Probes_Forward_A_transfer.genome_mfa.CT_conversion.fa.psl Filter.Probes_Forward_B_transfer.genome_mfa.CT_conversion_onlyCGprotection_reversecomplement.fa.psl Probes_Forward_A_transfer.genome_mfa.CT_conversion.fa.psl Probes_Forward_A_transfer.genome_mfa.CT_conversion_onlyCGprotection.fa.psl Probes_Forward_A_transfer.genome_mfa.CT_conversion_onlyCGprotection_reversecomplement.fa.psl Probes_Forward_A_transfer.genome_mfa.CT_conversion_reversecomplement.fa.psl Probes_Forward_B_transfer.genome_mfa.CT_conversion.fa.psl Probes_Forward_B_transfer.genome_mfa.CT_conversion_onlyCGprotection.fa.psl Probes_Forward_B_transfer.genome_mfa.CT_conversion_onlyCGprotection_reversecomplement.fa.psl Probes_Forward_B_transfer.genome_mfa.CT_conversion_reversecomplement.fa.psl

list=( $(ls | grep "^Probes_") )
for file in "${list[@]}";
 do rm Result/Filter.n.$file
 echo $file
 grep -P "\-\t\S+\t50\t\d+\t50\t\S+\t\d+\t\d+\t\d+\t1\t" $file > tem1.log
 grep -P "\-\t\S+\t51\t\d+\t51\t\S+\t\d+\t\d+\t\d+\t1\t" $file > tem2.log
 #grep -P "\-\t\S+\t50\t0\t\d+\t\S+\t\d+\t\d+\t\d+\t1\t" $file > tem1.log
 #grep -P "\-\t\S+\t51\t0\t\d+\t\S+\t\d+\t\d+\t\d+\t1\t" $file > tem2.log
 cat tem1.log tem2.log > tem.log
 rm tem1.log
 rm tem2.log
 while read -r line;
  do IFS='	' read -r -a array<<<"$line"
  #IFS=':' read -r -a arr<<<"${array[0]}"
  a1="${array[0]}";##match
  #echo $a1
  a11="${array[10]}";##qsize
  a2="${array[11]}";##qStart
  a3="${array[12]}";##qEnd
  identity=$(($a1*100/$a11)) #Browser window coordinates are 1-based [start,end] while PSL coordinates are 0-based [start,end)
  if [ "$a1" -ge "40" ]
   #then if [ "$identity" -ge "90" ];
    then echo "$line" >> Result/Filter.n.$file
    #echo $identity
   #fi
  fi
 done < "tem.log"
done
rm tem.log

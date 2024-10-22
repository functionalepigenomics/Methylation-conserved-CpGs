for filename in result_cis_tag-snpselection_in_*; do
sort $filename | awk '{duplicates[$1]++} END{for (ind in duplicates) {print ind,duplicates[ind]}}' > count_$filename
done 


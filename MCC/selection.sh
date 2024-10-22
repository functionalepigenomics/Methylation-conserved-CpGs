cat Filter.p.Probes_A.genome_mfa.GA_conversion.fa.psl Filter.n.Probes_A.genome_mfa.CT_conversion.fa.psl | cut -f10 | sort | uniq -u > a.log
cat Filter.n.Probes_B.genome_mfa.CT_conversion_onlyCGprotection.fa.psl Filter.p.Probes_B.genome_mfa.GA_conversion_onlyCGprotection.fa.psl | cut -f10 | sort | uniq -u > b.log
cat Filter.p.Probes_A.genome_mfa.CT_conversion_reversecomplement.fa.psl Filter.n.Probes_A.genome_mfa.GA_conversion_reversecomplement.fa.psl | cut -f10 | sort | uniq -u > c.log
comm -12 <( sort a.log ) <( sort b.log ) > d.log
comm -23 <( sort d.log ) <( sort c.log ) > probes_Gorilla_gorilla_gorilla.txt


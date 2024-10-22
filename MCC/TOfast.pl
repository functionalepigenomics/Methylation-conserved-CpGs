unlink("Probes_Forward_A.fasta");
unlink("Probes_Forward_B.fasta");
open FAM, "Probes_cut.txt";
open OUTA, ">>Probes_A.fasta";
open OUTB, ">>Probes_B.fasta";
while ($line=<FAM>){
chomp($line);
if ($line=~/^(\S+)\t(\S+)\t(\S+)/){
printf OUTA ">"."$1\n$2\n";
printf OUTB ">"."$1\n$3\n";
}
}

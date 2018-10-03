#!/bin/bash
# motif_analysis.sh
# 
# Seungsoo Kim
# October 3, 2018

module load perl/5.24.0
module load openmpi-1.10-x86_64
module load meme/4.12.0
module load bedtools/latest

sc="references/sacCer3.fa"
su="references/Sbay_revised.fa"

mcast --motif-pthresh 0.0005 --output-ethresh 2 --oc ../nobackup/motifs_cluster/Scer_0.0005_50_Leu3_Sdd4_Rgt1 motifs_Leu3_Sdd4_Rgt1.txt $sc
mcast --motif-pthresh 0.0005 --output-ethresh 2 --oc ../nobackup/motifs_cluster/Sbay_0.0005_50_Leu3_Sdd4_Rgt1 motifs_Leu3_Sdd4_Rgt1.txt $su

fimo --max-strand --thresh 0.0005 --o ../nobackup/motifs_Leu3/Scer_0.0005 motif_Leu3.txt $sc
fimo --max-strand --thresh 0.0005 --o ../nobackup/motifs_Sdd4/Scer_0.0005 motif_Sdd4.txt $sc
fimo --max-strand --thresh 0.0005 --o ../nobackup/motifs_Rgt1/Scer_0.0005 motif_Rgt1.txt $sc
fimo --max-strand --thresh 0.0005 --o ../nobackup/motifs_Leu3/Sbay_0.0005 motif_Leu3.txt $su
fimo --max-strand --thresh 0.0005 --o ../nobackup/motifs_Sdd4/Sbay_0.0005 motif_Sdd4.txt $su
fimo --max-strand --thresh 0.0005 --o ../nobackup/motifs_Rgt1/Sbay_0.0005 motif_Rgt1.txt $su

#bedtools intersect -c -a <(bedtools intersect -c -a <(bedtools intersect -c -a ../nobackup/motifs_cluster/Scer_0.0005_50_Leu3_Sdd4_Rgt1/mcast.gff -b ../nobackup/motifs_Leu3/Scer_0.0005/fimo.gff) -b ../nobackup/motifs_Sdd4/Scer_0.0005/fimo.gff) -b ../nobackup/motifs_Rgt1/Scer_0.0005/fimo.gff > Scer_0.0005_50_Leu3_Sdd4_Rgt1.txt
bedtools intersect -c -a <(bedtools intersect -c -a <(bedtools intersect -c -a ../nobackup/motifs_cluster/Scer_0.0005_50_Leu3_Sdd4_Rgt1/mcast.gff -b ../nobackup/motifs_Leu3/Scer_0.0005/fimo.gff) -b ../nobackup/motifs_Sdd4/Scer_0.0005/fimo.gff) -b ../nobackup/motifs_Rgt1/Scer_0.0005/fimo.gff | awk 'BEGIN{a = 1} $10 > 0 && $12 > 0 {print $0 "\t" a; a = a + 1}' > Scer_0.0005_50_Leu3_Sdd4_Rgt1_LR.txt
#bedtools intersect -c -a <(bedtools intersect -c -a <(bedtools intersect -c -a ../nobackup/motifs_cluster/Sbay_0.0005_50_Leu3_Sdd4_Rgt1/mcast.gff -b ../nobackup/motifs_Leu3/Sbay_0.0005/fimo.gff) -b ../nobackup/motifs_Sdd4/Sbay_0.0005/fimo.gff) -b ../nobackup/motifs_Rgt1/Sbay_0.0005/fimo.gff > Sbay_0.0005_50_Leu3_Sdd4_Rgt1.txt
bedtools intersect -c -a <(bedtools intersect -c -a <(bedtools intersect -c -a ../nobackup/motifs_cluster/Sbay_0.0005_50_Leu3_Sdd4_Rgt1/mcast.gff -b ../nobackup/motifs_Leu3/Sbay_0.0005/fimo.gff) -b ../nobackup/motifs_Sdd4/Sbay_0.0005/fimo.gff) -b ../nobackup/motifs_Rgt1/Sbay_0.0005/fimo.gff | awk 'BEGIN{a = 1} $10 > 0 && $12 > 0 {print $0 "\t" a; a = a + 1}' > Sbay_0.0005_50_Leu3_Sdd4_Rgt1_LR.txt

cat <(bedtools intersect -wa -wb -a Scer_0.0005_50_Leu3_Sdd4_Rgt1_LR.txt -b ../nobackup/motifs_Leu3/Scer_0.0005/fimo.gff) <(bedtools intersect -wa -wb -a Scer_0.0005_50_Leu3_Sdd4_Rgt1_LR.txt -b ../nobackup/motifs_Sdd4/Scer_0.0005/fimo.gff) <(bedtools intersect -wa -wb -a Scer_0.0005_50_Leu3_Sdd4_Rgt1_LR.txt -b ../nobackup/motifs_Rgt1/Scer_0.0005/fimo.gff) | awk -F"\t" '{OFS="\t"; print $1, $4, $5, $1 ":" $4 "-" $5, $6, $13, $17-$4, $18-$4, $19, $20, $22}' | sed 's/;qvalue.*//' | sed 's/Name.*Alias=//' | sed 's/;ID.*pvalue=/	/' > Scer_0.0005_50_Leu3_Sdd4_Rgt1_LR.motifs.txt
cat <(bedtools intersect -wa -wb -a Sbay_0.0005_50_Leu3_Sdd4_Rgt1_LR.txt -b ../nobackup/motifs_Leu3/Sbay_0.0005/fimo.gff) <(bedtools intersect -wa -wb -a Sbay_0.0005_50_Leu3_Sdd4_Rgt1_LR.txt -b ../nobackup/motifs_Sdd4/Sbay_0.0005/fimo.gff) <(bedtools intersect -wa -wb -a Sbay_0.0005_50_Leu3_Sdd4_Rgt1_LR.txt -b ../nobackup/motifs_Rgt1/Sbay_0.0005/fimo.gff) | awk -F"\t" '{OFS="\t"; print $1, $4, $5, $1 ":" $4 "-" $5, $6, $13, $17-$4, $18-$4, $19, $20, $22}' | sed 's/;qvalue.*//' | sed 's/Name.*Alias=//' | sed 's/;ID.*pvalue=/	/' > Sbay_0.0005_50_Leu3_Sdd4_Rgt1_LR.motifs.txt


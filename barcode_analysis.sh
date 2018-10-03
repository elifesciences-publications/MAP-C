#!/bin/bash
# barcode_analysis.sh
# Seungsoo Kim

# load modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load bowtie2/2.2.3 samtools/1.3

data='data'

name=$1
r1=$3
r2=$4
out='nobackup/'$5
annot=$6

# align reads using bowtie2
bowtie2 --very-sensitive -x nobackup/bowtie2/fixed-locus_ko_barcodes -U $data/$r1 2> $out/aligned/$name.1.out | samtools view -b - > $out/aligned/$name.1.bam

# count reads per reference barcode, and merge with annotations
join <(sort -k1,1 $annot) <(samtools view $out/aligned/$name.1.bam | awk '$5 >= 20 {OFS="\t"; print $3}' | sort | uniq -c | awk '{OFS="\t"; print $2, $1}') | awk '{for (i=2; i<NF; i++) printf i "\t"; print NF}' > $out/$name.counts.txt 


#!/bin/bash
# in-gene_ko.sh
# Seungsoo Kim

# load modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load bowtie2/2.2.3 samtools/1.3
module load python/2.7.3 cutadapt/1.8.3
module load bedtools/2.26.0

# fixed paths
primers='primers.txt'
data='data'
out='nobackup/in-gene_ko'
bt2='nobackup/bowtie2'

# sample info
name=$1
samptype=$2
r1=$3
r2=$4

# use cutadapt to trim adapter sequence
pr1='GATAGCTTCGGCGGCTAAGAC'
cutadapt -g $pr1 -o $out/cutadapt/$name.2.fastq.gz $data/$r2 > $out/cutadapt/$name.out

# use bowtie2 to map barcode sequences to expected barcodes
bowtie2 --very-sensitive -x $bt2/in-gene_ko_barcodes -U $data/$r1 2> $out/aligned/$name.1.out | samtools view -b - > $out/aligned/$name.1.bam

# tally count
paste <(samtools view $out/aligned/$name.1.bam | cut -f3,5) <(zcat $out/cutadapt/$name.2.fastq.gz | awk 'NR%4==2 {print substr($1,0,10)}') | awk '$2 > 0 {OFS="\t"; print $1, $3}' | sort -k1,1 -k2,2 | uniq -c | awk '{OFS="\t"; print $2, $3, $1}' | sort -k1,1 -k3,3nr -k2,2 > $out/$name.counts.withbcs.txt

awk 'BEGIN{OFS="\t"; a = ""; b = 0} {if (a == $1) b = b + $3; else {print a, b; a = $1; b = $3}} END{print a, b}' $out/$name.counts.withbcs.txt | tail -n +2 > $out/$name.counts.txt

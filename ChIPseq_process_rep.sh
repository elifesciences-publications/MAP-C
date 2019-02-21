#!/bin/bash
# ChIPseq_process_rep.sh
# Seungsoo Kim

# load modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load python/2.7.3 cutadapt/1.8.3
module load htslib/1.3.2 bowtie2/2.2.3 samtools/1.3.1_htslib-1.3.2
module load bedtools/2.26.0

# use cutadapt to trim adapters and low-quality ends of reads
out='nobackup/2019-01-20_ChIPseq/cutadapt'
in='data'
cutadapt -q 20 -m 28 -a AGATCGGAAGAGCACACGTC -A AGATCGGAAGAGCGTCGTGT -o $out/$1.1.fastq.gz -p $out/$1.2.fastq.gz $in/$2 $in/$3 > $out/$1.log

# use bowtie2 to align reads
in='nobackup/ChIPseq/cutadapt'
out='nobackup/ChIPseq/aligned'
bt2='nobackup/bowtie2'
bowtie2 --very-sensitive -X 2000 -p $NSLOTS -x $bt2/sacCer3 -1 $in/$1.1.fastq.gz -2 $in/$1.2.fastq.gz 2> $out/$1.log | samtools view -b - > $out/$1.bam

# insert sizes
samtools view $out/$1.bam | awk '$5 >= 30 && $9 > 0 && $9 <= 2000 {OFS="\t"; print $9}' | head -n 10000 > $out/$1.insertsizes.txt

# filter low mapping quality reads and pairs with >=1 read unmapped
samtools view -q 30 -F12 -@ $NSLOTS -b $out/$1.bam > $out/$1.filter.bam
samtools sort -@ $NSLOTS -o $out/$1.sorted.bam $out/$1.filter.bam
samtools rmdup $out/$1.sorted.bam $out/$1.dedup.bam

# calculate genome coverage using bedtools
in='nobackup/ChIPseq/aligned'
out='nobackup/ChIPseq/bedgraph'
genes='references/sacCer3_genes.gff'
genome='references/sacCer3.genome'
bedtools genomecov -ibam $in/$1.dedup.bam -g $genome -pc -bg -trackline > $out/$1.bedgraph
gzip $out/$1.bedgraph

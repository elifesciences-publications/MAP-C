#!/bin/bash
# in-gene_ko_all.sh
# Seungsoo Kim

# load modules
module load bowtie2/2.2.3 samtools/1.3
module load python/2.7.3 cutadapt/1.8.3
module load bedtools/2.26.0

# fixed paths
samples='in-gene_ko_samples.txt'
out='nobackup/in-gene_ko'
bt2='nobackup/bowtie2'

# create output directories
for x in cutadapt aligned sge
do
	mkdir -p $out/$x
done

mkdir -p $bt2

# index fasta files
samtools faidx in-gene_ko_barcodes.fa

# build bowtie2 indices
bowtie2-build in-gene_ko_barcodes.fa $bt2/in-gene_ko_barcodes

while read name samptype r1 r2
do
	echo $name
	qsub -N $name -o $out/sge -e $out/sge -l mfree=1G -cwd in-gene_ko.sh $name $samptype $r1 $r2
done < $samples

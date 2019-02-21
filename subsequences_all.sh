#!/bin/bash
# subsequences_all.sh
# Seungsoo Kim

# load modules
module load pear/0.9.5 bowtie2/2.2.3 samtools/1.3
module load python/2.7.3 cutadapt/1.8.3
module load bedtools/2.26.0

# fixed paths
samples='subsequences_samples.txt'
out='nobackup/subsequences'
bt2='nobackup/bowtie2'

# compile C++ code
g++ read_end_position.cpp -o read_end_position

# create output directories
for x in cutadapt cutadapt2 cutadapt3 paired aligned sge
do
	mkdir -p $out/$x
done

mkdir -p $bt2

# index fasta files
samtools faidx HAS1pr-TDA1pr.fa

# build bowtie2 indices
bowtie2-build HAS1pr-TDA1pr.fa $bt2/HAS1pr-TDA1pr

while read name samptype r1 r2
do
	mkdir -p $out/paired/$name
	echo $name
	qsub -N $name -o $out/sge -e $out/sge -l mfree=1G -pe serial 4 -cwd subsequences.sh $name $samptype $r1 $r2
done < $samples

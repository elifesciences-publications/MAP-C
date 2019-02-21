#!/bin/bash
# 3bp_substitution_all.sh
# Seungsoo Kim

# load modules
module load pear/0.9.5 bowtie2/2.2.3 samtools/1.3
module load python/2.7.3 cutadapt/1.8.3
module load bedtools/2.26.0

# fixed paths
samples='3bp_substitution_samples.txt'
out='nobackup/3bp_substitution'
bt2='nobackup/bowtie2'

# create output directories
for x in paired aligned sge
do
	mkdir -p $out/$x
done

mkdir -p $bt2

# index fasta files
samtools faidx 3bp_substitution.fa

# build bowtie2 indices
bowtie2-build 3bp_substitution.fa $bt2/3bp_substitution

while read name samptype r1 r2
do
	echo $name
#	qsub -q ravana.q -N $name -o $out/sge -e $out/sge -l mfree=1G -pe serial 4 -cwd 3bp_substitution.sh $name $ref $r1 $r2
	./3bp_substitution.sh $name $r1 $r2
done < $samples

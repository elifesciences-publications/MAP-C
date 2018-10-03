#!/bin/bash
# error-prone_all.sh
# Seungsoo Kim

# load modules
module load pear/0.9.5 bowtie2/2.2.3 samtools/1.3
module load python/2.7.3 cutadapt/1.8.3
module load bedtools/2.26.0

# fixed paths
samples='error-prone_samples.txt'
out='nobackup/error-prone'

# compile C++ code
g++ -std=c++11 annotate_mutations.cpp -o annotate_mutations

# create output directories
for x in paired aligned sge
do
	if [ ! -d $out/$x ]; then
		mkdir $out/$x
	fi
done

# index fasta files
samtools faidx error-prone.fa

# build bowtie2 indices
bowtie2-build error-prone.fa nobackup/bowtie2/error-prone

while read name samptype r1 r2
do
	echo $name
#	qsub -q ravana.q -N $name -o $out/sge -e $out/sge -l mfree=1G -pe serial 4 -cwd error-prone.sh $name $r1 $r2
	./error-prone.sh $name $r1 $r2
done < $samples

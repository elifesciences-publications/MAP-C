#!/bin/bash
# subsequences.sh
# Seungsoo Kim

# load modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load pear/0.9.5 bowtie2/2.2.3 samtools/1.3
module load python/2.7.3 cutadapt/1.8.3
module load bedtools/2.26.0

primers='subsequences_primers.txt'
data='data'
out='nobackup/subsequences'
bt2='nobackup/bowtie2'

# sample info
name=$1
samptype=$2
r1=$3
r2=$4

# calculate reverse complements of primer sequences
pr1=$(awk -v s=$samptype '$1==s {print $2}' $primers)
pr1rc=$(echo $pr1 | rev | tr 'ATCG' 'TAGC')
pr2=$(awk -v s=$samptype '$1==s {print $3}' $primers)
pr2rc=$(echo $pr2 | rev | tr 'ATCG' 'TAGC')
pr3=$(awk -v s=$samptype '$1==s {print $4}' $primers)
pr3rc=$(echo $pr3 | rev | tr 'ATCG' 'TAGC')

# trim primer sequences
cutadapt -u 4 -g $pr1 -G $pr2 -o $out/cutadapt/$name.1.fastq.gz -p $out/cutadapt/$name.2.fastq.gz $data/$r1 $data/$r2 > $out/cutadapt/$name.out
cutadapt -m 30 -a $pr2rc -A $pr1rc -o $out/cutadapt2/$name.1.fastq.gz -p $out/cutadapt2/$name.2.fastq.gz $out/cutadapt/$name.1.fastq.gz $out/cutadapt/$name.2.fastq.gz > $out/cutadapt2/$name.out

# use PEAR to merge read pairs
mkdir -p $out/paired/$name
if [ "$pr3" = "" ]; then
	gunzip $out/cutadapt2/$name.1.fastq.gz
	gunzip $out/cutadapt2/$name.2.fastq.gz
	pear -f $out/cutadapt2/$name.1.fastq -r $out/cutadapt2/$name.2.fastq -o $out/paired/$name > $out/paired/$name.out
	gzip $out/cutadapt2/$name.1.fastq
	gzip $out/cutadapt2/$name.2.fastq
else
	cutadapt -m 30 -A $pr3rc -g $pr3 -o $out/cutadapt3/$name.1.fastq.gz -p $out/cutadapt3/$name.2.fastq.gz $out/cutadapt2/$name.1.fastq.gz $out/cutadapt2/$name.2.fastq.gz > $out/cutadapt3/$name.out
	gunzip $out/cutadapt3/$name.1.fastq.gz
	gunzip $out/cutadapt3/$name.2.fastq.gz
	pear -f $out/cutadapt3/$name.1.fastq -r $out/cutadapt3/$name.2.fastq -o $out/paired/$name > $out/paired/$name.out
	gzip $out/cutadapt3/$name.1.fastq
	gzip $out/cutadapt3/$name.2.fastq	
fi

# align reads
bowtie2 --very-sensitive -x $bt2/HAS1pr-TDA1pr -U $out/paired/$name.assembled.fastq 2> $out/aligned/$name.out | samtools view -b - > $out/aligned/$name.bam

# extract read end positions to convert to bed file
samtools view $out/aligned/$name.bam | awk '$5 > 0 {OFS="\t"; print $3, $4, $6}' | ./read_end_position > $out/aligned/$name.bed

# use bedtools to calculate read coverage
bedtools genomecov -d -i $out/aligned/$name.bed -g HAS1pr-TDA1pr.genome > $out/aligned/$name.cov.txt


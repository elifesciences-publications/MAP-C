#!/bin/bash
# 3bp_substitution.sh
# Seungsoo Kim

# load modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load pear/0.9.5 bowtie2/2.2.3 samtools/1.3
module load python/2.7.3 cutadapt/1.8.3
module load bedtools/2.26.0

# fixed paths
samples='3bp_substitution_samples.txt'
data='data'
out='nobackup/3bp_substitution'

# sample info
name=$1
r1=$3
r2=$4
r1unzip=$(echo $r1 | sed 's/\.gz$//')
r2unzip=$(echo $r2 | sed 's/\.gz$//')

# use PEAR to merge reads and trim adapters
gunzip $data/$r1
gunzip $data/$r2
pear -j $NSLOTS -q 10 -f $data/$r1unzip -r $data/$r2unzip -o $out/paired/$name > $out/paired/$name.out
gzip $data/$r1unzip
gzip $data/$r2unzip
gzip $out/paired/$name.*.fastq

# create FASTA files from merged reads, with sequence names as rank_numreads
zcat $out/paired/$name.assembled.fastq.gz | awk 'NR%4==2 {print substr($1, 5, length($1)-4)}' | sort | uniq -c | awk '{OFS="\t"; print $2, $1}' | sort -k2,2nr | awk '{print ">" NR "_" $2; print $1}' > $out/paired/$name.counts.fa

# use Bowtie2 to align merged reads to reference
bowtie2 --very-sensitive -p $NSLOTS --reorder -f -x nobackup/bowtie2/$ref -U $out/paired/$name.counts.fa 2> $out/aligned/$name.out | samtools view -b - > $out/aligned/$name.bam

# annotate mutations - output contains 1) read count, 2) num mutations, 3) num substitutions, 4) num deletions, 5) num insertions, 6) mutation string
./annotate_mutations 3bp_substitution.fa <(samtools view $out/aligned/$name.bam) | sed 's/.*_//' | awk '{OFS="\t"; print $6, $2, $3, $4, $5, $1}' > $out/$name.mutcounts.txt

# count number of reads per substitution
python count_mutations.py <(awk '$3>=3 || $2==0 {OFS="\t"; print $1, $6; a = a + $6} END {print "All", a}' $out/$name.mutcounts.txt) | grep -v "D" | grep -v "+" > $out/$name.subcounts.txt


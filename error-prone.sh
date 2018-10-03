#!/bin/bash
# error-prone.sh
# Seungsoo Kim

# load modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load pear/0.9.5 bowtie2/2.2.3 samtools/1.3
module load python/2.7.3 cutadapt/1.8.3
module load bedtools/2.26.0
# fixed paths
samples='error-prone_samples.txt'
data='data'
out='nobackup/error-prone'
# sample info
name=$1
r1=$2
r2=$3
r1unzip=$(echo $r1 | sed 's/\.gz$//')
r2unzip=$(echo $r2 | sed 's/\.gz$//')

# use PEAR to merge paired reads
gunzip $data/$r1
gunzip $data/$r2
pear -j $NSLOTS -q 10 -f $data/$r1unzip -r $data/$r2unzip -o $out/paired/$name > $out/paired/$name.out
gzip $data/$r1
gzip $data/$r2
gzip $out/paired/$name.*.fastq

# create FASTA files from merged reads, with sequence names as rank_numreads
zcat $out/paired/$name.assembled.fastq.gz | awk 'NR%4==2 {print substr($1, 5, length($1)-4)}' | sort | uniq -c | awk '{OFS="\t"; print $2, $1}' | sort -k2,2nr | awk '{print ">" NR "_" $2; print $1}' > $out/paired/$name.counts.fa

# use Bowtie2 to align merged reads to reference
bowtie2 --very-sensitive -f -x nobackup/bowtie2/error-prone -U $out/paired/$name.counts.fa -S $out/aligned/$name.sam 2> $out/aligned/$name.out

# annotate mutations - output contains 1) read count, 2) num mutations, 3) num substitutions, 4) num deletions, 5) num insertions, 6) mutation string
./annotate_mutations error-prone.fa $out/aligned/$name.sam | sed 's/.*_//' | awk '{OFS="\t"; print $6, $2, $3, $4, $5, $1}' > $out/$name.mutcounts.txt

# read counts per mutation, filtering on various types of reads
python count_mutations.py <(awk '{OFS="\t"; print $1, $6; a = a + $6} END {print "All", a}' $out/$name.mutcounts.txt) | grep -v "D" | grep -v "+" > $out/$name.subcounts.txt


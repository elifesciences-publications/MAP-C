#!/bin/bash
# ScSuRNAseq_all.sh
# Seungsoo Kim

# load modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs

# fixed paths
out='nobackup/ScSuRNAseq'
bt2='nobackup/bowtie2'
ref='references/ScSu.fa'
samps='ScSuRNAseq_samples.txt'

# make output folders
for fold in 'cutadapt' 'aligned' 'counts' 'sge'
do
  mkdir -p $out'/'$fold
done

mkdir -p $bt2

# index fasta files
if [ ! -f $ref.fai ]; then
  samtools faidx $ref
fi

# build bowtie2 indices
if [ ! -f $bt2/ScSu.1.bt2 ]; then
  bowtie2-build $ref $bt2/ScSu
fi

# run samples through pipeline
while read samp read1 read2
do
    qsub -N $samp -o $out/sge -e $out/sge -l mfree=2G -pe serial 4 -cwd ./ScSuRNAseq.sh $samp $read1 $read2
done < $samps

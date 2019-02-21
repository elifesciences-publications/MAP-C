#!/bin/bash
# ChIPseq_process_rep_all.sh
# Seungsoo Kim

# load modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load bowtie2/2.2.3
module load bedtools/2.26.0

# fixed paths
out='nobackup/ChIPseq'
bt2='nobackup/bowtie2'
ref='references/sacCer3.fa'
samps='ChIPseq_samples.txt'

# make output folders
for fold in 'sge' 'cutadapt' 'aligned' 'bedgraph'
do
  mkdir -p $out'/'$fold
done

mkdir -p $bt2

# make bowtie2 reference
if [ ! -f $bt2/sacCer3.1.bt2 ]; then
  bowtie2-build $ref $bt2/sacCer3
fi

# run samples through pipeline
while read samp read1 read2
do
    qsub -N $samp -o $out/sge -e $out/sge -l mfree=2G -pe serial 2 -cwd ./ChIPseq_process_rep.sh $samp $read1 $read2
done < $samps

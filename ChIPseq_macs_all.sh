#!/bin/bash
# ChIPseq_macs_all
# Seungsoo Kim

out='nobackup/ChIPseq'

for strain in 'Leu3' 'Sdd4' 'Rgt1'
do
	for cond in 'saturated' 'exponential'
	do
		qsub -N ${strain}_${cond}_merged -o $out/sge -e $out/sge -l mfree=2G -cwd ./ChIPseq_macs.sh ${strain}_${cond}_IP_1.dedup.bam ${strain}_${cond}_IP_2.dedup.bam ${strain}_${cond}_IP_3.dedup.bam ${strain}_${cond}_input_1.dedup.bam ${strain}_${cond}_input_2.dedup.bam ${strain}_${cond}_input_3.dedup.bam ${strain}_${cond}_merged
	done
done

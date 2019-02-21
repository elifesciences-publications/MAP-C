#!/bin/bash
# ChIPseq_annotate_genes.sh 
# Seungsoo Kim

# load modules
module load bedtools/2.26.0

peaksdir='nobackup/ChIPseq/macs'

for tf in 'Leu3' 'Sdd4' 'Rgt1'
do
	for c in 'saturated' 'exponential'
	do
		bedtools closest -D a -fd -t first -a <(bedtools closest -D a -fu -t first -a <(awk '{OFS="\t"; print $1, $2+$10, $2+$10+1, $7}' $peaksdir/${tf}_${c}_merged_peaks.narrowPeak | sort -k1,1 -k2,2n) -b <(sort -k1,1 -k4,4n ../sacCer3_genes_cleaned.gff | cut -f1,4-10)) -b <(sort -k1,1 -k4,4n ../sacCer3_genes_cleaned.gff | cut -f1,4-10) | cut -f1-4,9,11-13,18,20-22 > $peaksdir/${tf}_${c}_merged_peaks_annotated.bed

		awk '{OFS="\t"; if ($5=="-" && $9=="-") print $0, $7; else if ($5=="+" && $9=="+") print $0, $11; else if ($5=="-" && $9=="+" && (-$8 < $12)) print $0, $7; else if ($5=="-" && $9=="+" && (-$8 >= $12)) print $0, $11; else print $0, $7 "-" $11}' $peaksdir/${tf}_${c}_merged_peaks_annotated.bed > $peaksdir/${tf}_${c}_merged_peaks_labeled.bed
	done
done

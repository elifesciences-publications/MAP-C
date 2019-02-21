#!/bin/bash
# ChIPseq_compare_conditions.sh
# Seungsoo Kim

# load modules
module load perl/5.24.0
module load openmpi-1.10-x86_64
module load meme/4.12.0
module load bedtools/2.26.0

peakdir='nobackup/ChIPseq/macs'

# count motifs within 100 bp of summit
for tf in 'Leu3' 'Sdd4' 'Rgt1'
do
	bedtools intersect -loj -wa -wb -a <(awk '{OFS="\t"; print $1, $2+$10-100, $2+$10+101, $7}' $peakdir/${tf}_saturated_merged_peaks.narrowPeak) -b <(awk '{OFS="\t"; print $1, $2+$10-100, $2+$10+101, $7}' $peakdir/${tf}_exponential_merged_peaks.narrowPeak) | cut -f1-4,8 | awk '{OFS="\t"; if ($5 == ".") print $1, $2, $3, $4, "0"; else print $0}' > $peakdir/${tf}_saturated_merged_overlap_exponential_peaks.bed
	bedtools intersect -loj -wa -wb -a <(awk '{OFS="\t"; print $1, $2+$10-100, $2+$10+101, $7}' $peakdir/${tf}_exponential_merged_peaks.narrowPeak) -b <(awk '{OFS="\t"; print $1, $2+$10-100, $2+$10+101, $7}' $peakdir/${tf}_saturated_merged_peaks.narrowPeak) | cut -f1-4,8 | awk '{OFS="\t"; if ($5 == ".") print $1, $2, $3, $4, "0"; else print $0}' > $peakdir/${tf}_exponential_merged_overlap_saturated_peaks.bed
done


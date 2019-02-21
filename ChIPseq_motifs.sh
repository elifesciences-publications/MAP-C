#!/bin/bash
# ChIPseq_motifs.sh
# Seungsoo Kim

# load modules
module load perl/5.24.0
module load openmpi-1.10-x86_64
module load meme/4.12.0
module load bedtools/2.26.0

# paths
peakdir='nobackup/ChIPseq/macs'
out='nobackup/ChIPseq/meme'
mkdir -p $out

for tf in 'Leu3' 'Sdd4' 'Rgt1'
do
	for c in 'saturated' 'exponential'
	do
		# MEME for de novo motif discovery
		bedtools getfasta -fi references/sacCer3.fa -bed <(awk '{OFS="\t"; print $1, $2+$10-50, $2+$10+50}' $peakdir/${tf}_${c}_merged_peaks.narrowPeak) -fo $peakdir/${tf}_${c}_w50.fa
		bedtools getfasta -fi references/sacCer3.fa -bed <(awk '$7 >= 2 {OFS="\t"; print $1, $2+$10-50, $2+$10+50}' $peakdir/${tf}_${c}_merged_peaks.narrowPeak) -fo $peakdir/${tf}_${c}_gt2_w50.fa
		bedtools getfasta -fi references/sacCer3.fa -bed <(sort -k7,7nr $peakdir/${tf}_${c}_merged_peaks.narrowPeak | head -n 50 | awk '{OFS="\t"; print $1, $2+$10-50, $2+$10+50}') -fo $peakdir/${tf}_${c}_top50_w50.fa
		bedtools getfasta -fi references/sacCer3.fa -bed <(sort -k7,7nr $peakdir/${tf}_${c}_merged_peaks.narrowPeak | head -n 100 | awk '{OFS="\t"; print $1, $2+$10-50, $2+$10+50}') -fo $peakdir/${tf}_${c}_top100_w50.fa
		meme -dna -revcomp -nmotifs 3 -minw 6 -maxw 10 -oc $out/${tf}_${c}_top100_w50 $peakdir/${tf}_${c}_top100_w50.fa
		meme -dna -revcomp -nmotifs 3 -minw 6 -maxw 10 -oc $out/${tf}_${c}_top50_w50 $peakdir/${tf}_${c}_top50_w50.fa
		meme -dna -revcomp -nmotifs 3 -minw 6 -maxw 10 -oc $out/${tf}_${c}_w50 $peakdir/${tf}_${c}_w50.fa
		meme -dna -revcomp -nmotifs 3 -minw 6 -maxw 10 -oc $out/${tf}_${c}_gt2_w50 $peakdir/${tf}_${c}_gt2_w50.fa
		meme -dna -revcomp -nmotifs 3 -minw 6 -maxw 20 -oc $out/${tf}_${c}_top100_w50_max20 $peakdir/${tf}_${c}_top100_w50.fa
		meme -dna -revcomp -nmotifs 3 -minw 6 -maxw 20 -oc $out/${tf}_${c}_top50_w50_max20 $peakdir/${tf}_${c}_top50_w50.fa
		meme -dna -revcomp -nmotifs 3 -minw 6 -maxw 20 -oc $out/${tf}_${c}_w50_max20 $peakdir/${tf}_${c}_w50.fa
		meme -dna -revcomp -nmotifs 3 -minw 6 -maxw 20 -oc $out/${tf}_${c}_gt2_w50_max20 $peakdir/${tf}_${c}_gt2_w50.fa

		# count motifs within 100 bp of summit
		p='0.001'
		bedtools window -w 100 -c -a <(awk '{OFS="\t"; print $1, $2+$10, $2+$10+1, $7}' $peakdir/${tf}_${c}_merged_peaks.narrowPeak) -b <(cut -f1,4,5,6,7,8 nobackup/motifs_${tf}/Scer_${p}/fimo.gff) > $peakdir/${tf}_${c}_merged_summits_motifs_$p.bed
		bedtools window -w 100 -c -a <(bedtools window -w 100 -c -a <(bedtools window -w 100 -c -a <(awk '{OFS="\t"; print $1, $2+$10, $2+$10+1, $7}' $peakdir/${tf}_${c}_merged_peaks.narrowPeak) -b <(cut -f1,4,5,6,7,8 nobackup/motifs_Leu3/Scer_${p}/fimo.gff)) -b <(cut -f1,4-8 nobackup/motifs_Sdd4/Scer_${p}/fimo.gff)) -b <(cut -f1,4-8 nobackup/motifs_Rgt1/Scer_${p}/fimo.gff) > $peakdir/${tf}_${c}_merged_summits_motifs_all_$p.bed
	done
done

# extract tRNA gene sequences, +/- 50 bp from center
bedtools getfasta -fi references/sacCer3.fa -bed <(awk '{OFS="\t"; print $1, int(($4+$5)/2)-50, int(($4+$5)/2)+50}' sacCer3_tRNA.gff) -fo sacCer3_tRNA_w50.fa
# MEME analysis of tRNA genes
meme -dna -revcomp -nmotifs 3 -minw 6 -maxw 20 -oc $out/sacCer3_tRNA_cent sacCer3_tRNA_w50.fa

# annotate with gene names
for tf in 'Leu3' 'Sdd4' 'Rgt1'
do
    for c in 'saturated' 'exponential'
    do  
        bedtools closest -D a -fd -t first -a <(bedtools closest -D a -fu -t first -a <(awk '{OFS="\t"; print $1, $2+$10, $2+$10+1, $7}' $peaksdir/${tf}_${c}_merged_peaks.narrowPeak | sort -k1,1 -k2,2n) -b <(sort -k1,1 -k4,4n sacCer3_genes_cleaned.gff | cut -f1,4-10)) -b <(sort -k1,1 -k4,4n sacCer3_genes_cleaned.gff | cut -f1,4-10) | cut -f1-4,9,11-13,18,20-22 > $peaksdir/${tf}_${c}_merged_peaks_annotated.bed

        awk '{OFS="\t"; if ($5=="-" && $9=="-") print $0, $7; else if ($5=="+" && $9=="+") print $0, $11; else if ($5=="-" && $9=="+" && (-$8 < $12)) print $0, $7; else if ($5=="-" && $9=="+" && (-$8 >= $12)) print $0, $11; else print $0, $7 "-" $11}' $peaksdir/${tf}_${c}_merged_peaks_annotated.bed > $peaksdir/${tf}_${c}_merged_peaks_labeled.bed
    done
done

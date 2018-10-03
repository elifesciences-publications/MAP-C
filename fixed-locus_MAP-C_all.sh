#!/bin/bash
# fixed-locus_MAP-C_all.sh
# Seungsoo Kim
# September 28, 2018

# load modules
module load bowtie2/2.2.3 samtools/1.3
module load python/2.7.3 cutadapt/1.8.3
module load bedtools/2.26.0

# index fasta files
samtools faidx fixed-locus_ko_barcodes.fa

# build bowtie2 indices
bowtie2-build fixed-locus_ko_barcodes.fa nobackup/bowtie2/fixed-locus_ko_barcodes

# fixed paths
for expt in 'fixed-locus_ko' 'validation_ko' 'interactor_ko' 'Rgt1_10aaD'
do
	out='nobackup/'$expt
	samples=${expt}_samples.txt

	# create output directories
	for x in aligned sge
	do
		mkdir -p $out/$x
	done

	while read name samptype r1 r2
	do
		echo $name
		qsub -q ravana.q -N $name -o $out/sge -e $out/sge -l mfree=1G -cwd barcode_analysis.sh $name $samptype $r1 $r2 $expt ${expt}_annotations.txt
	done < $samples
done

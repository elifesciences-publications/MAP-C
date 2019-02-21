#!/bin/bash
# motif_analysis.sh
# Seungsoo Kim

# load modules
module load perl/5.24.0
module load openmpi-1.10-x86_64
module load meme/4.12.0
module load bedtools/2.26.0

# paths
mcast='nobackup/mcast'
macs='nobackup/ChIPseq/macs'
mkdir -p $mcast
sc='references/sacCer3.fa'
su='references/Sbay_revised.fa'

# process S. cer genes
grep -v Dubious references/sacCer3_genes.gff | awk '{OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8, $9}' | sed 's/ID=//' | sed 's/;Name=.*gene=/  /' | sed 's/;Name=.*//' | sed 's/;Alias.*//' | awk '{OFS="\t"; if (NF==10) print $0; else print $0, $9}' > references/sacCer3_genes_cleaned.gff

# process S. bay (uvarum) genes
join -1 9 <(cat references/Sbay_revised_genes.gff | sed 's/ID=.*BLAST=//' | sed 's/;ncbi.*//' | awk '{OFS="\t"; if ((substr($9,1,1) == "Y") && (length($9) <= 9)) print $0}' | sort -k9,9) <(cut -f9,10 ../sacCer3_genes_cleaned.gff | sort -k1,1) | awk '{OFS="\t"; print $2, $3, $4, $5, $6, $7, $8, $9, $1, $10}' | sort -k1,1V -k4,4n > references/Sbay_revised_genes_cleaned.gff

# run MCAST
for p in '0.001' '0.0005'
do
	for g in '50' '100' '200'
	do
		mcast --motif-pthresh $p --max-gap $g --output-ethresh 2 --oc $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1 motifs_Leu3_Sdd4_Rgt1.txt $sc
		mcast --motif-pthresh $p --max-gap $g --output-ethresh 2 --oc $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1 motifs_Leu3_Sdd4_Rgt1.txt $su
	done
done

# run FIMO for each TF
for tf in 'Leu3' 'Sdd4' 'Rgt1'
do
	mkdir -p nobackup/motifs_$tf
	for p in '0.001' '0.0005'
	do
		fimo --max-strand --thresh $p --oc nobackup/motifs_$tf/Scer_$p motif_$tf.txt $sc
		fimo --max-strand --thresh $p --oc nobackup/motifs_$tf/Sbay_$p motif_$tf.txt $su
	done
done

# overlap MCAST hits with FIMO hits
for p in '0.001' '0.0005'
do
	for g in '50' '100' '200'
	do
		bedtools intersect -c -a <(bedtools intersect -c -a <(bedtools intersect -c -a $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1/mcast.gff -b nobackup/motifs_Leu3/Scer_${p}/fimo.gff) -b nobackup/motifs_Sdd4/Scer_${p}/fimo.gff) -b nobackup/motifs_Rgt1/Scer_${p}/fimo.gff > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1.txt
		bedtools intersect -c -a <(bedtools intersect -c -a <(bedtools intersect -c -a $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1/mcast.gff -b nobackup/motifs_Leu3/Sbay_${p}/fimo.gff) -b nobackup/motifs_Sdd4/Sbay_${p}/fimo.gff) -b nobackup/motifs_Rgt1/Sbay_${p}/fimo.gff > $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1.txt
	done
done

# annotate left and right side genes
for p in '0.001' '0.0005'
do
    for g in '50' '100' '200'
    do  
        bedtools closest -D a -fd -t first -a <(bedtools closest -D a -fu -t first -a <(sort -k1,1 -k4,4n $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1.txt | cut -f1,4,5,6-12) -b <(sort -k1,1 -k4,4n references/sacCer3_genes_cleaned.gff | cut -f1,4-10)) -b <(sort -k1,1 -k4,4n references/sacCer3_genes_cleaned.gff | cut -f1,4-10) | cut -f1-4,7-10,15,17-19,24,26-28 | sort -k4,4nr > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_annotated.txt
        bedtools closest -D a -fd -t first -a <(bedtools closest -D a -fu -t first -a <(sort -k1,1 -k4,4n $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1.txt | cut -f1,4,5,6-12) -b <(sort -k1,1 -k4,4n references/Sbay_revised_genes_cleaned.gff | cut -f1,4-10)) -b <(sort -k1,1 -k4,4n references/Sbay_revised_genes_cleaned.gff | cut -f1,4-10) | cut -f1-4,7-10,15,17-19,24,26-28 | sort -k4,4nr > $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_annotated.txt

	    # extract only homologous, in SGD gene name order
		# homologous on L side
		join -1 10 -2 10 <(sort -k10,10 $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_annotated.txt) <(sort -k10,10 $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_annotated.txt) | awk '{OFS="\t"; print $2, $3, $4, $5, $6, $7, $8, $9, $10, $1, $11, $12, $13, $14, $15, $16}' > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologousL.txt
		# homologous on R side
		join -1 14 -2 14 <(sort -k14,14 $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_annotated.txt) <(sort -k14,14 $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_annotated.txt) | awk '{OFS="\t"; print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $1, $15, $16}' > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologousR.txt
		# merge
		cat $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologousL.txt $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologousR.txt | sort -k10,10 | uniq > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous.txt

		# homologous on L side
		join -1 10 -2 10 <(sort -k10,10 $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_annotated.txt) <(sort -k10,10 $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_annotated.txt) | awk '{OFS="\t"; print $2, $3, $4, $5, $6, $7, $8, $9, $10, $1, $11, $12, $13, $14, $15, $16}' > $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologousL.txt
		# homologous on R side
		join -1 14 -2 14 <(sort -k14,14 $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_annotated.txt) <(sort -k14,14 $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_annotated.txt) | awk '{OFS="\t"; print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $1, $15, $16}' > $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologousR.txt
		# merge
		cat $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologousL.txt $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologousR.txt | sort -k10,10 | uniq > $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous.txt

		# label clusters as promoters
		awk '{OFS="\t"; if ($9=="-" && $13=="-") print $0, $11 "pr"; else if ($9=="+" && $13=="+") print $0, $15 "pr"; else if ($9=="-" && $13=="+" && (-$12 < $16)) print $0, $11 "pr"; else if ($9=="-" && $13=="+" && (-$12 >= $16)) print $0, $15 "pr"; else print $0, $11 "-" $15}' $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous.txt > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_labeled.txt
		awk '{OFS="\t"; if ($9=="-" && $13=="-") print $0, $11 "pr"; else if ($9=="+" && $13=="+") print $0, $15 "pr"; else if ($9=="-" && $13=="+" && (-$12 < $16)) print $0, $11 "pr"; else if ($9=="-" && $13=="+" && (-$12 >= $16)) print $0, $15 "pr"; else print $0, $11 "-" $15}' $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous.txt > $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_labeled.txt

		awk '{OFS="\t"; if ($9=="-" && $13=="-") print $0, $11 "pr"; else if ($9=="+" && $13=="+") print $0, $15 "pr"; else if ($9=="-" && $13=="+" && (-$12 < $16)) print $0, $11 "pr"; else if ($9=="-" && $13=="+" && (-$12 >= $16)) print $0, $15 "pr"; else print $0, $11 "-" $15}' $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_annotated.txt > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_labeled.txt
		awk '{OFS="\t"; if ($9=="-" && $13=="-") print $0, $11 "pr"; else if ($9=="+" && $13=="+") print $0, $15 "pr"; else if ($9=="-" && $13=="+" && (-$12 < $16)) print $0, $11 "pr"; else if ($9=="-" && $13=="+" && (-$12 >= $16)) print $0, $15 "pr"; else print $0, $11 "-" $15}' $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_annotated.txt > $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_labeled.txt

		# number clusters, individually
		awk '{OFS="\t"; print $0, NR}' $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_labeled.txt > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_numbered.txt
		awk '{OFS="\t"; print $0, NR}' $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_labeled.txt > $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_numbered.txt

		# number clusters, only homologous
		join -1 17 -2 17 <(sort -k17,17 $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_labeled.txt) <(sort -k17,17 $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_labeled.txt) | awk '{print $0, $5 + $21}' | sort -k34,34nr | awk '{OFS="\t"; print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $1, NR}' > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_numbered.txt
		join -1 17 -2 17 <(sort -k17,17 $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_labeled.txt) <(sort -k17,17 $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_labeled.txt) | awk '{print $0, $5 + $21}' | sort -k34,34nr | awk '{OFS="\t"; print $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $1, NR}' > $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_numbered.txt
		join -1 17 -2 17 <(sort -k17,17 $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_labeled.txt) <(sort -k17,17 $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_labeled.txt) | awk '{print $0, $5 + $21}' | sort -k34,34nr | awk '$9 > 0 && $7 > 0 && $8 > 0' | awk '{OFS="\t"; print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $1, NR}' > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR_numbered.txt
		join -1 17 -2 17 <(sort -k17,17 $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_labeled.txt) <(sort -k17,17 $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_labeled.txt) | awk '{print $0, $5 + $21}' | sort -k34,34nr | awk '$9 > 0 && $7 > 0 && $8 > 0' | awk '{OFS="\t"; print $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $1, NR}' > $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR_numbered.txt

		# extract FIMO motif hits overlapping motif clusters
		cat <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_numbered.txt -b nobackup/motifs_Leu3/Scer_${p}/fimo.gff) <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_numbered.txt -b nobackup/motifs_Sdd4/Scer_${p}/fimo.gff) <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_numbered.txt -b nobackup/motifs_Rgt1/Scer_${p}/fimo.gff) | awk -F"\t" '{OFS="\t"; print $1, $2, $3, $1 ":" $2 "-" $3, $4, $17, $18, $22-$2, $23-$2, $24, $25, $27}' | sed 's/;qvalue.*//' | sed 's/Name.*Alias=//' | sed 's/;ID.*pvalue=/	/' > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous.motifs.txt
		cat <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_numbered.txt -b nobackup/motifs_Leu3/Sbay_${p}/fimo.gff) <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_numbered.txt -b nobackup/motifs_Sdd4/Sbay_${p}/fimo.gff) <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_numbered.txt -b nobackup/motifs_Rgt1/Sbay_${p}/fimo.gff) | awk -F"\t" '{OFS="\t"; print $1, $2, $3, $1 ":" $2 "-" $3, $4, $17, $18, $22-$2, $23-$2, $24, $25, $27}' | sed 's/;qvalue.*//' | sed 's/Name.*Alias=//' | sed 's/;ID.*pvalue=/	/' > $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous.motifs.txt

		cat <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR_numbered.txt -b nobackup/motifs_Leu3/Scer_${p}/fimo.gff) <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR_numbered.txt -b nobackup/motifs_Sdd4/Scer_${p}/fimo.gff) <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR_numbered.txt -b nobackup/motifs_Rgt1/Scer_${p}/fimo.gff) | awk -F"\t" '{OFS="\t"; print $1, $2, $3, $1 ":" $2 "-" $3, $4, $17, $18, $22-$2, $23-$2, $24, $25, $27}' | sed 's/;qvalue.*//' | sed 's/Name.*Alias=//' | sed 's/;ID.*pvalue=/	/' > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR.motifs.txt
		cat <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR_numbered.txt -b nobackup/motifs_Leu3/Sbay_${p}/fimo.gff) <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR_numbered.txt -b nobackup/motifs_Sdd4/Sbay_${p}/fimo.gff) <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR_numbered.txt -b nobackup/motifs_Rgt1/Sbay_${p}/fimo.gff) | awk -F"\t" '{OFS="\t"; print $1, $2, $3, $1 ":" $2 "-" $3, $4, $17, $18, $22-$2, $23-$2, $24, $25, $27}' | sed 's/;qvalue.*//' | sed 's/Name.*Alias=//' | sed 's/;ID.*pvalue=/	/' > $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR.motifs.txt

		cat <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_numbered.txt -b nobackup/motifs_Leu3/Scer_${p}/fimo.gff) <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_numbered.txt -b nobackup/motifs_Sdd4/Scer_${p}/fimo.gff) <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_numbered.txt -b nobackup/motifs_Rgt1/Scer_${p}/fimo.gff) | awk -F"\t" '{OFS="\t"; print $1, $2, $3, $1 ":" $2 "-" $3, $4, $17, $18, $22-$2, $23-$2, $24, $25, $27}' | sed 's/;qvalue.*//' | sed 's/Name.*Alias=//' | sed 's/;ID.*pvalue=/	/' > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1.motifs.txt
		cat <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_numbered.txt -b nobackup/motifs_Leu3/Sbay_${p}/fimo.gff) <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_numbered.txt -b nobackup/motifs_Sdd4/Sbay_${p}/fimo.gff) <(bedtools intersect -wa -wb -F 0.6 -a $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1_numbered.txt -b nobackup/motifs_Rgt1/Sbay_${p}/fimo.gff) | awk -F"\t" '{OFS="\t"; print $1, $2, $3, $1 ":" $2 "-" $3, $4, $17, $18, $22-$2, $23-$2, $24, $25, $27}' | sed 's/;qvalue.*//' | sed 's/Name.*Alias=//' | sed 's/;ID.*pvalue=/	/' > $mcast/Sbay_${p}_${g}_Leu3_Sdd4_Rgt1.motifs.txt

		# extract maximum ChIP-seq signal (fold enrichment over input) for each TF
		rm $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_maxFE.txt
		rm $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR_maxFE.txt
		rm $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_maxFE.txt
		for s in 'Leu3' 'Sdd4' 'Rgt1'
		do
			for c in 'exponential' 'saturated'
			do
				bedtools intersect -wa -wb -a $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_numbered.txt -b $macs/${s}_${c}_merged_FE.bdg.gz | python get_max.py 18 22 | cut -f1-18,22 > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_${s}_${c}_maxFE.txt
				cut -f17-19 $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_${s}_${c}_maxFE.txt | awk -v s=$s -v c=$c '{OFS="\t"; print $0, s, c}' >> $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_maxFE.txt
				bedtools intersect -wa -wb -a $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR_numbered.txt -b $macs/${s}_${c}_merged_FE.bdg.gz | python get_max.py 18 22 | cut -f1-18,22 > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR_${s}_${c}_maxFE.txt
				cut -f17-19 $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR_${s}_${c}_maxFE.txt | awk -v s=$s -v c=$c '{OFS="\t"; print $0, s, c}' >> $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_homologous_LSR_maxFE.txt
				bedtools intersect -wa -wb -a $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_numbered.txt -b $macs/${s}_${c}_merged_FE.bdg.gz | python get_max.py 18 22 | cut -f1-18,22 > $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_${s}_${c}_maxFE.txt
				cut -f17-19 $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_${s}_${c}_maxFE.txt | awk -v s=$s -v c=$c '{OFS="\t"; print $0, s, c}' >> $mcast/Scer_${p}_${g}_Leu3_Sdd4_Rgt1_maxFE.txt
			done
		done
    done
done

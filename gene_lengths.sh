#!/bin/bash
# gene_lengths.sh
# Seungsoo Kim

dir='data'
out='nobackup'
mkdir -p $out

sed 's/;.*//' $dir/sacCer3_genes.gff | sed 's/ID=//' | awk '{OFS="\t"; print $9, $5-$4}' > $out/sacCer3_genes_length.txt
sed 's/;Parent=.*//' data/ScSu_genes.gff | sed 's/ID=.*;Gene=//' | sed 's/;Name=.*//' | sed 's/Gene=//' | awk '{OFS="\t"; print $9, $5-$4}' > $out/ScSu_genes_length.txt


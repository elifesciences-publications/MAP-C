#!/bin/bash
# ChIPseq_macs.sh
# Seungsoo Kim

# load modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load python/2.7.3 numpy/1.8.1 setuptools/25.1.1 MACS/2.1.0

# paths
bams='nobackup/ChIPseq/aligned'
bgs='nobackup/ChIPseq/bedgraph'
out='nobackup/ChIPseq/macs'

# call peaks
macs2 callpeak -t $bams/$1 $bams/$2 $bams/$3 -c $bams/$4 $bams/$5 $bams/$6 -f BAMPE -g 11860000 -n $7 --outdir $out/ -q 0.01 -B --SPMR --keep-dup all

# get fold enrichment track
macs2 bdgcmp -t $out/$7_treat_pileup.bdg -c $out/$7_control_lambda.bdg -o $out/$7_FE.bdg -m FE


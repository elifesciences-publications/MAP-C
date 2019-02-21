#!/bin/bash
# yetfasco.sh
# Seungsoo Kim

# path with PFMs
pfms='yetfasco'

# file with list of PFM names
pfmlist='yetfasco.txt'

while read x
do
	python motif_rc.py $pfms/$x.pfm > $pfms/${x}_rc.pfm
	python motif_pad.py $pfms/$x.pfm 15 > $pfms/${x}_pad15.pfm
	python motif_pad.py $pfms/${x}_rc.pfm 15 > $pfms/${x}_rc_pad15.pfm
done < $pfmlist

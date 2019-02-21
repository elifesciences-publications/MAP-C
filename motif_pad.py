# motif_pad.py
# Generates PFM padded up to n columns with equal frequency of A, C, T, G.
# USAGE: python motif_pad.py FILE N_COLUMNS
# 
# Seungsoo Kim
# February 20, 2019

import sys

f = open(sys.argv[1])
cols = int(sys.argv[2]) + 1

lines = f.readlines()
ncols = len(lines[0].strip().split('\t'))
for line in lines:
	if (cols >= ncols):
		print line.strip() + "\t0.25" * (cols-ncols)
	else:
		print line.strip()

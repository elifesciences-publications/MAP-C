# motif_rc.py
# Generates reverse complemented PFM.
# USAGE: python motif_rc.py FILE
# 
# Seungsoo Kim
# February 20, 2019

import sys

f = open(sys.argv[1])
lines = f.readlines()
for line in lines:
	line = line.strip().split('\t')
	base = line[0]
	del line[0]
	if base == 'A':
		print 'T' + "\t" + "\t".join(line[::-1])
	elif base == 'T':
		print 'A' + "\t" + "\t".join(line[::-1])
	elif base == 'C':
		print 'G' + "\t" + "\t".join(line[::-1])
	elif base == 'G':
		print 'C' + "\t" + "\t".join(line[::-1])


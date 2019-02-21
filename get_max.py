# get_max.py
# Given a tab-delimited stdin with at least 2 columns,
# extracts the maximum value of VAR_COL for each value 
# of FIX_COL
# USAGE: python get_max.py FIX_COL VAR_COL
#
# Seungsoo Kim
# February 20, 2019

import sys

#f=open(sys.argv[1])
fixcol=int(sys.argv[1])
varcol=int(sys.argv[2])
curfix=0
curvar=0
curline=""
lineno=0
for line in sys.stdin:
	line=line.strip()
	fields=line.split("\t")
	if (fields[fixcol-1]==curfix):
		if (fields[varcol-1] > curvar):
			curfix=fields[fixcol-1]
			curvar=fields[varcol-1]
			curline=line
	else:
		if (lineno != 0):
			print curline
		curfix=fields[fixcol-1]
		curvar=fields[varcol-1]
		curline=line
	lineno += 1
print curline

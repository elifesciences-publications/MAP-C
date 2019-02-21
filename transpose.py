# transpose.py
# Transposes a tab-delimited file, turning columns into
# rows and vice versa.
# USAGE: python transpose.py FILE
# 
# Seungsoo Kim
# February 20, 2019

import sys

f = open(sys.argv[1],'rb')
parsed = []
for line in f:
	filtered = filter(None, line.strip().split('\t'))
	parsed.append(filtered)

for j in range(len(parsed[0])):
	for i in range(len(parsed)):
		print parsed[i][j],
	print

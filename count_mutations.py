# count_mutations.py
#
# Seungsoo Kim
# September 27, 2018

import sys

mutfile = open(sys.argv[1])
mutdict = {}
for line in mutfile:
	line = line.strip().split("\t")
	muts = line[0].split(",")
	count = int(line[1])
	for mut in muts:
		if mut in mutdict:
			mutdict[mut] += count
		else:
			mutdict[mut] = count
for key in sorted(mutdict):
	print key + "\t" + str(mutdict[key])

import matplotlib
matplotlib.use('pdf')
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

parser = argparse.ArgumentParser(prog='plot-distances.py')
parser.add_argument('distances', metavar="DISTANCES", help="comma-separated list of files produced by bedtools closest output (-d).")
parser.add_argument('outfile', metavar="OUTFILE", help="name of the output file.")
args = parser.parse_args()

histograms = []
labels = []
infiles = args.distances.split(',')

for file in infiles:
	histogram = []
	if 'shuffled' in file:
		labels.append('shuffled assembly variants')
	else:
		labels.append('assembly variants')
	for line in open(file, 'r'):
		if line.startswith('#'):
			continue
		fields = line.split()
		histogram.append(int(fields[-1]))
	histograms.append(histogram)

plt.xlabel('distance to closest gnomAD variant (bp)')
plt.ylabel('Count')
plt.yscale('symlog')
plt.xscale('symlog')
#plt.xlim(0,10**10)
bins = np.logspace(-1, 8, 60)
bins = np.insert(bins, 0, 0.0, axis=0)
plt.hist(histograms, label=labels, bins=bins)
plt.legend(loc='upper right')
plt.savefig(args.outfile)

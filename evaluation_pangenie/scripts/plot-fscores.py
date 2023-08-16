import sys
import argparse
from collections import defaultdict
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt 
import math

import glob

class Stats:
	def __init__(self, filename):
		self._fscore = 0
		for line in open(filename, 'r'):
			if not "None" in line:
				continue
			fields = line.strip().split()
			assert len(fields) == 8
			assert fields[0] == "None"
			self._fscore = float(fields[7])
			
	def get_fscore(self):
		return self._fscore


def plot_fscores(mode, run_mode, variants, folder, outfile, coverage, sample, regions):
	methods = ['pangenie', 'pangenie'] if regions not in ['repeats', 'nonrep', 'external'] else ['pangenie', 'pangenie', 'bayestyper', 'graphtyper']
	legend_methods = ['PanGenie (all)', 'PanGenie (high-gq)'] if regions not in ['repeats', 'nonrep', 'external'] else ['PanGenie (all)', 'PanGenie (high-gq)', 'BayesTyper', 'GraphTyper']

	qualities = ['0', '200'] if regions not in ['repeats', 'nonrep', 'external'] else ['0', '200', '0', '0', '0', '0', '0']
#	assert len(coverages) < 6
#	sorted_coverages = sorted([int(coverage) for coverage in coverages])
#	coverages = [str(c) for c in sorted_coverages]

	markers = ['D', 'x', 'P', '*', 'v']
	plot_index = 1
	colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#808080']
	assert len(variants) <= 12
	n_rows = min(math.ceil(len(variants) / 3), 4)
	n_rows = 4
	n_cols = min(3, len(variants))
	n_cols = 3
	plt.figure(figsize=(13,20))
	plot_lines = []
	for varianttype in variants:#
		if "complex" in varianttype:
			assert not 'external' in regions
			region = regions + '-complex'
			var = varianttype.split('-')[0]
		else:
			if not 'external' in regions:
				region = regions + '-simple'
			else:
				region = regions
			var = varianttype
		if any(['discovery' in r for r in run_mode]):
			assert varianttype.startswith('sv') or varianttype.startswith('indel')
		labels = []
		considered_methods = []
		considered_colors = []
		fscore_all = []
		patterns = []
		for i, method, quality in zip(range(len(methods)), methods, qualities):
			# paragraph types only midsize and large variants. Graphtyper in all variants
			if (method in ['paragraph']) and (varianttype.startswith('snp') or varianttype.startswith('small') or varianttype.startswith('indel')):
				continue
			# gatk does not type large variants
			if (method == 'gatk') and (varianttype.startswith('large') or varianttype.startswith('sv')):
				continue
#			if any(['discovery' in r for r in run_mode]) and not method in ['gatk', 'platypus']:
#				continue
			for run in run_mode:
				if 'discovery' in run and not method in ['gatk', 'platypus']:
					continue
				considered_methods.append(legend_methods[i] + ' (discovery)' if 'discovery' in run else legend_methods[i])
				patterns.append("/" if 'discovery' in run else "")
				considered_colors.append(colors[i])
				filename= folder + sample + "/evaluation/precision-recall-typable/" + mode + "/" + method + "-" + run + '/coverage-' + str(coverage) + "_" + region + "_" + var + "/qual_" + quality + "/summary.txt"
				stats = Stats(filename)
				fscore = stats.get_fscore()
				fscore_all.append(fscore)
#		fig = plt.figure()
		ax = plt.subplot(n_rows, n_cols, plot_index)
		for m,f,c,p in zip(considered_methods, fscore_all, considered_colors,patterns):
			ax.bar(m,f,color=c,edgecolor='black', hatch=p)
#		plt.bar(considered_methods, fscore_all,color=considered_colors, hatch=patterns)
		plt.title(varianttype)
		plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.3)
		plt.ylabel('adjusted F-score')
		plt.xticks(rotation=45, ha='right')
		plt.ylim(0.0, 1.0)
		plot_index += 1
	plt.subplots_adjust(top=0.96, bottom=0.08, left=0.10, right=0.90, hspace=0.2, wspace=0.3)
	if sample is not None:
		plt.suptitle('Results for ' + sample, fontsize=20)
	plt.tight_layout()
	plt.savefig(outfile)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='plot-precision-recall.py', description=__doc__)
	parser.add_argument('mode', metavar='MODE', help='truth set (external|graph).')
	parser.add_argument('regions', metavar='REGIONS', help='regions to evaluate (repeats|non-repeats|external|kir|mhc).')
	parser.add_argument('-run_mode', metavar='RUNMODE', nargs='+', help='(discovery-biallelic).')
	parser.add_argument('-coverages', metavar='COVERAGE', help='coverages.')
#	parser.add_argument('-variants', metavar='VARIANTS', required=True, nargs='+', help='comma separated list of variants.')
	parser.add_argument('-folder', metavar='FOLDER', default='', help='path to folder with evaluation results per method.')
	parser.add_argument('-outfile', metavar='OUTPUT', default='qualities.pdf', help='name of the output-file')
	parser.add_argument('-sample', metavar='SAMPLE', required=True, help='name of the sample (used for title).')
	
	args = parser.parse_args()
	variants = []
	if args.mode == "external":
		if args.regions == "external":
			variants = ["indel", "sv"]
		else:
			variants = ["indel", "large-deletion", "large-insertion", "large-complex"]
	elif any(['discovery' in i for i in args.run_mode]):
		variants = ["indel", "indel-complex"]
	elif args.regions  not in ['repeats', 'nonrep', 'external']:
		variants = ["all", "all-complex"]
	else:
		variants = ["small-deletion", "small-insertion", "small-complex", "midsize-deletion", "midsize-insertion", "midsize-complex", "large-deletion", "large-insertion", "large-complex", "snp", "snp-complex"]
	coverage = args.coverages
	plot_fscores(args.mode, args.run_mode, variants, args.folder, args.outfile, coverage, args.sample, args.regions)

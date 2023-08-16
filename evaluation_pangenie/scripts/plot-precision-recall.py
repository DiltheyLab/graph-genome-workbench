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
		self._prec = None
		self._rec = None
		for line in open(filename, 'r'):
			if not "None" in line:
				continue
			fields = line.strip().split()
			assert len(fields) == 8
			assert fields[0] == "None"
			self._prec = float(fields[5])
			self._rec = float(fields[6])
			
	def get_precision(self):
		return self._prec
	def get_recall(self):
		return self._rec


def plot_prec_rec(mode, run_mode, variants, folder, outfile, coverages, sample, regions, variantset):
	methods = ['pangenie', 'pangenie'] if regions not in ['repeats', 'nonrep', 'external'] else ['pangenie', 'pangenie', 'bayestyper', 'graphtyper'] 
	legend_methods = ['PanGenie (all)', 'PanGenie (high-gq)'] if regions not in ['repeats', 'nonrep', 'external'] else ['PanGenie (all)', 'PanGenie (high-gq)', 'BayesTyper', 'GraphTyper']

	qualities = ['0', '200'] if regions not in ['repeats', 'nonrep', 'external'] else ['0', '200', '0', '0', '0', '0', '0']
	assert len(coverages) < 6
	#sorted_coverages = sorted([int(coverage) for coverage in coverages])
	coverages = [str(c) for c in coverages]

	markers = ['D', 'x', 'P', '*', 'v']
	plot_index = 1
	colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#808080']
	assert len(variants) <= 12
	n_rows = 4 #min(math.ceil(len(variants) / 3), 4)
	n_cols = 3 #min(3, len(variants))
	plt.figure(figsize=(13,20))
	plot_lines = []
	considered_methods = []
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
		if any(['discovery' in m for m in run_mode]):
			assert varianttype.startswith('sv') or varianttype.startswith('indel')
		plt.subplot(n_rows, n_cols, plot_index)
		labels = []
		for i, method, quality in zip(range(len(methods)), methods, qualities):
			# paragraph types only midsize and large variants. Avoid graphtyper
			if (method in ['paragraph']) and (varianttype.startswith('snp') or varianttype.startswith('small') or varianttype.startswith('indel')):
				continue
			# gatk does not type large variants
			if (method == 'gatk') and (varianttype.startswith('large') or varianttype.startswith('sv')):
				continue
#			if any(['discovery' in m for m in run_mode]) and not method in ['gatk', 'platypus']:
#				continue
			considered_methods.append(method)
			for run in run_mode:
				if 'discovery' in run and not method in ['gatk', 'platypus']:
					continue
				precision_all = []
				recall_all = []
				linestyle= '--' if (len(run_mode)>1) and ('discovery' in run) else '-'
				for j,coverage in enumerate(coverages):
					filename= folder + sample + "/evaluation/precision-recall-" + variantset + "/" + mode + "/" + method + "-" + run + '/coverage-' + str(coverage) + "_" + region + "_" + var + "/qual_" + quality + "/summary.txt"
					stats = Stats(filename)
					unique_kmers = []
					precision = stats.get_precision()
					recall = stats.get_recall()
					precision_all.append(precision)
					recall_all.append(recall)
					plt.scatter([recall], [precision], marker='o', s=float(j*13)+13, color=colors[i], zorder=2)
				label=legend_methods[i]
				plt.plot(recall_all, precision_all, linestyle=linestyle, label=label, color=colors[i], linewidth=0.8, zorder=1)
		plt.title(varianttype)
		plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.3)
		assert variantset in ['all', 'typable']
		x_label = 'recall' if variantset == 'all' else 'adjusted recall' 
		y_label = 'precision' if variantset == 'all' else 'adjusted precision'
		plt.xlabel(x_label)
		plt.ylabel(y_label)
#		plt.ylim([40,106])
#		plt.xlim([0,106])
		plot_index += 1
	plt.subplots_adjust(top=0.96, bottom=0.08, left=0.10, right=0.90, hspace=0.2, wspace=0.3)
	if sample is not None:
		plt.suptitle('Results for ' + sample, fontsize=20)
	# create legend
	handles = []
	for i, method in enumerate(methods):
		for run in run_mode:
			if not method in considered_methods:
				continue
			label = (legend_methods[i] + ' (discovery)') if 'discovery' in run else legend_methods[i]
			linestyle = '--' if (len(run_mode)>1) and ('discovery' in run) else '-'
			line = matplotlib.lines.Line2D([],[], color=colors[i], markersize=100, linewidth=2.0, linestyle=linestyle, label=label)
			handles.append(line)
	for i,coverage in enumerate(coverages):
		label = str(coverage) + 'x'
		line = matplotlib.lines.Line2D([],[], color='grey', marker='o', linestyle='none', label=label)
		handles.append(line)
	labels = [h.get_label() for h in handles]
	plt.legend(handles=handles, labels=labels)
	plt.savefig(outfile)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='plot-precision-recall.py', description=__doc__)
	parser.add_argument('mode', metavar='MODE', help='truth set (external|graph).')
	parser.add_argument('regions', metavar='REGIONS', help='regions to evaluate (repeats|non-repeats|external|kir|mhc).')
	parser.add_argument('-run_mode', metavar='RUNMODE', nargs='+', help='(genotyping-biallelic|discovery-biallelic).')
	parser.add_argument('-coverages', metavar='COVERAGES', nargs='+', help='coverages.')
#	parser.add_argument('-variants', metavar='VARIANTS', required=True, nargs='+', help='comma separated list of variants.')
	parser.add_argument('-folder', metavar='FOLDER', default='', help='path to folder with evaluation results per method.')
	parser.add_argument('-outfile', metavar='OUTPUT', default='qualities.pdf', help='name of the output-file')
	parser.add_argument('-sample', metavar='SAMPLE', required=True, help='name of the sample (used for title).')
	parser.add_argument('-variantset', metavar="VARIANTSET", default="typable", help="(all|typable)")
	
	args = parser.parse_args()
	variants = []
	if args.mode == "external":
		if args.regions == 'external':
			variants = ["indel", "sv"]
		else:
			variants = ["indel", "large-deletion", "large-insertion", "large-complex"]
	elif any(['discovery' in i for i in args.run_mode]):
		variants = ["indel", "indel-complex"]
	elif args.regions  not in ['repeats', 'nonrep', 'external']:
		variants = ["all", "all-complex"]
	else:
		variants = ["small-deletion", "small-insertion", "small-complex", "midsize-deletion", "midsize-insertion", "midsize-complex", "large-deletion", "large-insertion", "large-complex", "snp", "snp-complex"]
	coverages = args.coverages
	plot_prec_rec(args.mode, args.run_mode, variants, args.folder, args.outfile, coverages, args.sample, args.regions, args.variantset)

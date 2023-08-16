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
		self._conc = None
		self._typed = None
		for line in open(filename, 'r'):
			fields = line.strip().split()
			if fields[0] == 'weighted_concordance':
				assert len(fields) == 2
				assert fields[0] == 'weighted_concordance'
				self._conc = float(fields[1]) * 100.0
				continue
			if fields[0] == 'typed':
				assert len(fields) == 2
				assert fields[0] == 'typed'
				self._typed = float(fields[1]) * 100.0
	def get_concordance(self):
		return self._conc
	def get_typed(self):
		return self._typed


def plot_concordance(mode, run_mode, variants, folder, outfile, coverages, sample, regions):
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
	print(n_rows,n_cols)
	plt.figure(figsize=(13,20))
	plot_lines = []
	considered_methods = []
	for varianttype in variants:
		if "complex" in varianttype:
			assert regions != "external"
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
			# paragraph types only midsize and large variants
			if (method in ['paragraph']) and (varianttype.startswith('snp') or varianttype.startswith('small') or varianttype.startswith('indel')):
				continue
			# graphtyper types only SVs. Graphtyper in all--> comment all
			#if (method == "graphtyper") and (varianttype.startswith('snp') or varianttype.startswith('small') or varianttype.startswith('midsize') or varianttype.startswith('indel')):
				#continue
			# gatk does not type large variants
			if (method == 'gatk') and (varianttype.startswith('large') or varianttype.startswith('sv')):
				continue
#			if any(['discovery' in m for m in run_mode]) and not method in ['gatk', 'platypus']:
#				continue
			considered_methods.append(method)
			for run in run_mode:
				if 'discovery' in run and not method in ['gatk', 'platypus']:
					continue
				concordance_all = []
				typed_all = []
				linestyle= '--' if (len(run_mode)>1) and ('discovery' in run) else '-'
				for j,coverage in enumerate(coverages):
					filename= folder + sample + "/evaluation/concordance/" + mode + "/" + method + "-" + run + '/coverage-' + str(coverage) + "_" + region + "_" + var + "/qual_" + quality + "/summary.txt"
					stats = Stats(filename)
					unique_kmers = []
					concordance = stats.get_concordance()
					typed = stats.get_typed()
					concordance_all.append(concordance)
					typed_all.append(typed)
					plt.scatter([typed], [concordance], marker='o', color=colors[i], s=float(j*13)+13, zorder=2)
				label=(legend_methods[i] + ' (discovery)') if 'discovery' in run else legend_methods[i]
				plt.plot(typed_all, concordance_all, linestyle=linestyle, label=label, color=colors[i], linewidth=0.8, zorder=1)
		plt.title(varianttype)
		plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.3)
		plt.xlabel('genotyped variants [%]')
		plt.ylabel('weighted genotype concordance [%]')
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
	parser = argparse.ArgumentParser(prog='plot-concordances.py', description=__doc__)
	parser.add_argument('mode', metavar='MODE', help='truth set (external|graph).')
	parser.add_argument('regions', metavar='REGIONS', help='regions to evaluate (repeats|non-repeats|external|kir|mhc).')
	parser.add_argument('-coverages', metavar='COVERAGES', nargs='+', help='coverages.')
	parser.add_argument('-run_mode', metavar='RUNMODE', nargs='+',  help='(genotyping-biallelic|discovery-biallelic|both).')
#	parser.add_argument('-variants', metavar='VARIANTS', required=True, nargs='+', help='comma separated list of variants.')
	parser.add_argument('-folder', metavar='FOLDER', default='', help='path to folder with evaluation results per method.')
	parser.add_argument('-outfile', metavar='OUTPUT', default='qualities.pdf', help='name of the output-file')
	parser.add_argument('-sample', metavar='SAMPLE', required=True, help='name of the sample (used for title).')
	
	args = parser.parse_args()
	coverages = args.coverages
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
	plot_concordance(args.mode, args.run_mode, variants, args.folder, args.outfile, coverages, args.sample, args.regions)

import sys
import argparse
import math

def parse_concordance(filename):
	concordance = None
	typed = None
	for line in open(filename, 'r'):
		if line.startswith('weighted_concordance'):
			concordance = float(line.split()[-1]) * 100.0
		if line.startswith('typed'):
			typed = float(line.split()[-1]) * 100.0
	assert concordance is not None
	assert typed is not None
	return concordance, typed


def parse_precision_recall(filename):
	precision = None
	recall = None
	fscore = None
	for line in open(filename, 'r'):
		if 'None' in line:
			fields = line.strip().split()
			assert len(fields) == 8
			assert fields[0] == 'None'
			precision = float(fields[5]) if fields[5] != 'None'  else 0.0
			recall = float(fields[6]) if fields[6] != 'None' else 0.0
			fscore = float(fields[7]) if fields[7] != 'None' else 0.0
			for k in [precision, recall, fscore]:
				k = 0.0 if math.isnan(k) else k
	assert precision is not None
	assert recall is not None
	assert fscore is not None
	return precision, recall, fscore


def collect_concordances(folder, variant, outname, coverage, sample, regions, filt, version):
	with open(outname + '.tsv', 'w') as outfile:
		outfile.write('\t'.join(['sample', 'weighted_genotype_concordance', 'typed_variants']) + '\n')
		filename = folder + '/' + version + '-' + coverage + '-' + filt + '-' + regions + '-' + variant + '/concordance/summary.txt'
		concordance, typed = parse_concordance(filename)
		outfile.write('\t'.join([sample, str(concordance), str(typed)]) + '\n')


def collect_precision_recall(folder, variant, outname, coverage, sample, regions, filt, version):
	with open(outname + '.tsv', 'w') as outfile:
		outfile.write('\t'.join(['sample', 'precision', 'recall', 'fscore']) + '\n')
		filename = folder + '/' + version + '-' + coverage + '-' + filt + '-' + regions + '-' + variant + '/precision-recall/summary.txt'
		precision, recall, fscore = parse_precision_recall(filename)
		outfile.write('\t'.join([sample, str(precision), str(recall), str(fscore)]) + '\n')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='collect-results.py', description=__doc__)
	parser.add_argument('metric', metavar='METRIC', help='concordance|precision-recall')
	parser.add_argument('version', metavar='VERSIONS', help='genotyping algorithm used')
	parser.add_argument('sample', metavar='SAMPLE', help='comma-separated list of samples.')
	parser.add_argument('coverage', metavar='COVERAGE', help='coverage')
	parser.add_argument('region', metavar='REGIONS', help='region to evaluate')
	parser.add_argument('filter', metavar='FILTERS', help='filter to evaluate')
	parser.add_argument('-variants', metavar='VARIANTS', required=True, help='variant type.')
	parser.add_argument('-folder', metavar='FOLDER', default='', help='path to folder with evaluation results per method.')
	parser.add_argument('-outfile', metavar='OUTPUT', default='qualities.pdf', help='name of the output-file')

	args = parser.parse_args()
	assert args.metric in ['concordance', 'precision-recall']
	sample = args.sample
	if args.metric == 'concordance':
		collect_concordances(args.folder, args.variants, args.outfile, args.coverage, sample, args.region, args.filter, args.version)
	if args.metric == 'precision-recall':
		collect_precision_recall(args.folder, args.variants, args.outfile, args.coverage, sample, args.region, args.filter, args.version)

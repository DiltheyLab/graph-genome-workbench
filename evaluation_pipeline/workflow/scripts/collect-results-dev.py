import sys
import argparse

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


def parse_precision_recall_vcfeval(filename):
	precision = None
	recall = None
	fscore = None
	for line in open(filename, 'r'):
		if 'None' in line:
			fields = line.strip().split()
			assert len(fields) == 8
			assert fields[0] == 'None'
			precision = float(fields[5])
			recall = float(fields[6])
			fscore = float(fields[7])
	assert precision is not None
	assert recall is not None
	assert fscore is not None
	return precision, recall, fscore

def parse_precision_recall_truvari(filename):
	precision = None
	recall = None
	fscore = None
	with open(filename) as json_file:
		data = json.load(json_file)
	
	precision = data['precision'] if data['precision'] is not None else 0
	recall = data['recall'] if data['recall'] is not None else 0
	f1score = data['f1'] if data['f1'] is not None else 0
	assert precision is not None
	assert recall is not None
	assert fscore is not None
	return precision, recall, fscore


def collect_concordances(folder, variant, outname, coverage, samples, regions):
	with open(outname + '.tsv', 'w') as outfile:
		outfile.write('\t'.join(['sample', 'weighted_genotype_concordance', 'typed_variants']) + '\n')
		for i, sample in enumerate(samples):
			filename = folder + '/' + sample + '/' + coverage + '/concordance/' + regions + '_' + variant + '/summary.txt'
			concordance, typed = parse_concordance(filename)
			outfile.write('\t'.join([sample, str(concordance), str(typed)]) + '\n')


def collect_precision_recall(folder, variant, outname, coverage, samples, regions):
	with open(outname + '.tsv', 'w') as outfile:
		outfile.write('\t'.join(['sample', 'precision', 'recall', 'fscore']) + '\n')
		for i, sample in enumerate(samples):
			filename = folder + '/' + sample + '/' + coverage + '/precision-recall-typable/' + regions + '_' + variant + '/summary.txt'
			if metric == 'precision-recall-vcfeval':
				precision, recall, fscore = parse_precision_recall_vcfeval(filename)
			elif metric == 'precision-recall-truvari':
				precision, recall, fscore = parse_precision_recall_truvari(filename)
			else:
				assert(False)
			outfile.write('\t'.join([sample, str(precision), str(recall), str(fscore)]) + '\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='collect-results.py', description=__doc__)
	parser.add_argument('metric', metavar='METRIC', help='concordance|precision-recall-vcfeval|precision-recall-truvari')
	parser.add_argument('coverage', metavar='COVERAGE', help='coverage')
	parser.add_argument('sample', metavar='SAMPLE', help='Sample')
	parser.add_argument('regions', metavar='REGIONS', help='region to evaluate')
	parser.add_argument('-variants', metavar='VARIANTS', required=True, help='variant type.')
	parser.add_argument('-folder', metavar='FOLDER', default='', help='path to folder with evaluation results per method.')
	parser.add_argument('-outfile', metavar='OUTPUT', default='qualities.pdf', help='name of the output-file')

	args = parser.parse_args()
	assert args.metric in ['concordance', 'precision-recall-vcfeval', 'precision-recall-truvari']
	samples = [args.sample]
	if args.metric == 'concordance':
		collect_concordances(args.folder, args.variants, args.outfile, args.coverage, samples, args.regions)
	else:
		collect_precision_recall(args.folder, args.variants, args.outfile, args.coverage, samples, args.regions)
	
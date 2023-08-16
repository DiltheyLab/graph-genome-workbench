import sys
import argparse
import matplotlib.pyplot as plt

def parse_log(filename):
	names = filename.split('/')
	assert len(names) >= 9 
	sample = names[-8]
	assert names[-3].startswith('coverage-')
	region = names[-3].split('_')[-2]

	baseline = 0
	ids = 0
	for line in open(filename, 'r'):
		if 'baseline' in line:
			baseline = int(line.split()[1])
		if 'ids' in line:
			ids = int(line.split()[1])
	return sample, region, baseline, ids

def parse_file(filename):
	# determine region this file corresponds to
	# coverage-30_hladpa1-simple_all/qual_0/summary.txt
	names = filename.split('/')
	assert len(names) >= 9 
	sample = names[-8]	
	assert names[-3].startswith('coverage-')
	region = names[-3].split('_')[-2]


	number_absent = 0
	number_het = 0
	number_hom = 0
	conc_absent = 0.0
	conc_het = 0.0
	conc_hom = 0.0
	conc_weight = 0.0
	
	# read concordances
	for line in open(filename, 'r'):
		if line.startswith('0') and (len(line.split()) == 5):
			number_absent = sum([int(i) for i in line.split()[1:]])
		if line.startswith('1'):
			number_het = sum([int(i) for i in line.split()[1:]])
		if line.startswith('2'):
			number_hom = sum([int(i) for i in line.split()[1:]])
		if line.startswith('genotype_concordance_absent'):
			conc_absent = float(line.split()[-1])
		if line.startswith('genotype_concordance_het'):
			conc_het = float(line.split()[-1])
		if line.startswith('genotype_concordance_hom'):
			conc_hom = float(line.split()[-1])
		if line.startswith('weighted_concordance'):
			conc_weight = float(line.split()[-1])
#	print(number_absent, number_het, number_hom)
#	print(conc_absent, conc_het, conc_hom)
	non_zero = 0
	conc = 0
	for c,n in zip([conc_absent, conc_het, conc_hom], [number_absent, number_het, number_hom]):
		if not n > 0:
			continue
		non_zero += 1
		conc += c
	final_conc = conc / non_zero if non_zero > 0 else 0
	if non_zero == 3:
		assert abs(final_conc - conc_weight) <= 0.0000000001
	return sample, region, final_conc

parser = argparse.ArgumentParser(prog='plot-concordances.py', description=__doc__)
parser.add_argument('files', metavar='FILES', nargs='+',  help='genotype concordance files')
parser.add_argument('--outname', metavar='OUTNAME', default='plot.pdf', help='name of output file')
args = parser.parse_args()

samples = set([])
regions = set([])

sample_to_numbers = {}
sample_to_stats = {}

for filename in args.files:
	if filename.endswith('.txt'):
		sample, region, conc = parse_file(filename)
		samples.add(sample)
		regions.add(region)
		sample_to_numbers[(sample,region)] = conc
	else:
		sample,region,baseline,ids = parse_log(filename)
		sample_to_stats[(sample, region)] = [baseline, ids]

for region in sorted(regions):
	sum = 0
	for sample in ['HG00731', 'NA24385', 'NA12878']: #sorted(samples):
		print('\t'.join([sample, region, str(sample_to_numbers[(sample,region)]*100.0), str(sample_to_stats[(sample,region)][1]) + ' (' + str(sample_to_stats[(sample,region)][0])  + ')']))
		sum += sample_to_numbers[(sample,region)]
	print('average\t\t' + str((sum/len(samples))*100.0))
	sample_to_numbers[('average', region)] = sum/len(samples)


# plotting
regions = sorted(regions)
x_values = []
i=0
num=0
while num < len(regions):
	i += 1
	if i % 3 == 0:
		continue
	x_values.append(i)
	num += 1

region_to_name = {
	'hlaa' : 'HLA-A',
	'hlab' : 'HLA-B',
	'hlac' : 'HLA-C',
	'hladpa1' : 'HLA-DPA1',
	'hladpb1' : 'HLA-DPB1',
	'hladra' : 'HLA-DRA',
	'hladrb1' : 'HLA-DRB1',
	'hladqa1' : 'HLA-DQA1',
	'hladqb1' : 'HLA-DQB1',
	'c4': 'C4',
	'hladrb5': 'HLA-DRB5'
}

print(x_values)
bubble_to_name = {'complex': 'complex', 'simple':'biallelic'}
markers = ['o', 'v', 'P']
regions_labels = [region_to_name[r.split('-')[0]] + '-' + bubble_to_name[r.split('-')[1]] for r in regions]
print(regions_labels)
plt.grid(True)
for i,sample in enumerate(list(samples) + ['average']):
	if sample == 'average':
		marker = '_'
		alpha = 1.0
		color='red'
	else:
		marker = markers[i]
		alpha = 1.0
		color='black'
	plt.scatter(x_values, [sample_to_numbers[(sample,i)]*100.0 for i in regions], alpha=alpha, color=color, marker=marker, label=sample)
plt.xticks(x_values, regions_labels, rotation=90)
plt.legend(loc="lower left")
plt.ylabel('weighted genotype concordance [%]')
plt.tight_layout()
plt.savefig(args.outname)

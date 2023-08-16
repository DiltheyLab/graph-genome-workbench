#!/usr/bin/python

import argparse
import math
import sys
import random
import re
from collections import defaultdict
from variantclassifier import VariantType, determine_variant_type

def select_samples(samples_to_choose, number):
	selected_sample_ids=random.sample(samples_to_choose, number)
	selected_sample_ids.sort(key=lambda tup: tup[1])
	return selected_sample_ids


parser = argparse.ArgumentParser(prog='sample-panels.py', description=__doc__)
parser.add_argument('sample', metavar='sample', help='sample to genotype (will be excluded from panel).')
parser.add_argument('numbers', metavar='nr_samples', help='comma separated list of numbers of samples to choose.', type=str)
parser.add_argument('outname', metavar='outname', help='prefix of name of output reference panel VCF file.')
parser.add_argument('--truth', metavar='truth', help='also generate ground truth VCF containing genotypes of given sample and write it to this file.', default=None)
args = parser.parse_args()

total_variants = defaultdict(int)
skipped_variants = defaultdict(int)

numbers=[int(i) for i in args.numbers.split(',')]
truth_file=None
if args.truth:
	# ground truth output file
	truth_file = open(args.truth, 'w')
# output files
panel_filenames=[args.outname + '.' + str(n) + '.vcf' for n in numbers]
panel_files=[open(filename, 'w') for filename in panel_filenames]
# pool of selected samples to generate VCFs from
selected_sample_ids=None
# id of truth sample
sample_id=None
# last position seen
prev_chrom=None
prev_pos=0
# statisitcs
skip_not_typable=0
skip_absent=0
skip_overlapping=0
skip_missing=0

for line in sys.stdin:
	if line.startswith('##'):
		# description line
		if args.truth:
			truth_file.write(line)
		for file in panel_files:
			file.write(line)
		continue
	if line.startswith('#'):
		# header line
		if selected_sample_ids != None:
			raise ValueError('Malformatted input VCF file.')
		fields=line.split()
		# all vcf samples
		all_samples = fields[9:]
		# randomly choose args.number samples
		samples_to_choose=[]
		for i,s in enumerate(all_samples):
			if s == args.sample:
				sample_id = i + 9
			else:
				samples_to_choose.append( (s,i + 9) )
		if len(samples_to_choose) < max(numbers):
			raise ValueError('Not enough samples given.')
		selected_sample_ids=select_samples(samples_to_choose, max(numbers))
		sys.stderr.write('Selected samples: ' + ','.join([s[0] for s in selected_sample_ids]) + '\n')
		# output header line
		header_line=fields[:9]
		for file,number in zip(panel_files,numbers):
			file.write('\t'.join(header_line + [s[0] for s in selected_sample_ids[:number]]) + '\n')
		if args.truth:
			truth_file.write('\t'.join(header_line) + '\t' + args.sample + '\n')
		continue
	# vcf line
	if selected_sample_ids == None:
		raise ValueError('Malformatted input VCF file.')
	if sample_id == None:
		raise ValueError('Given sample not contained in VCF.')
	fields=line.split()
	skip_line=False
	skip_truth=False
	chrom = fields[0]
	position = fields[1]
	# get IDs from INFO field
	variant_type = determine_variant_type(line)
	total_variants[variant_type] += 1
	# skip implicit variants
	if any(c not in 'CAGTcagt,' for c in fields[4]) or any(c not in 'CAGTcagt,' for c in fields[3]):
		sys.stderr.write('Skipping position ' + chrom + ':' + position + ' since alternative allele is not given explicitly.\n')
		continue
	# set filter to PASS
	fields[6] = 'PASS'
	for file,number in zip(panel_files,numbers):
		vcf_line = fields[:9]
		for s in selected_sample_ids[:number]:
			vcf_line.append(fields[s[1]])
		alleles = '|'.join(vcf_line[9:9+min(numbers)]).split('|')
		true_alleles = fields[sample_id].split('|')
		# determine how many alleles are missing
		if '.' in true_alleles:
			skip_missing += 1
			skip_line = True
			skip_truth = True
			continue
		# skip positions absent in all paths (wouldn't be called as variants given only the reference paths)
		if all(a == '0' or a == '.' for a in alleles):
			sys.stderr.write('Skipping position ' + chrom + ':' + position + ' since variant is absent in all samples.\n')
			skip_absent += 1
			skip_line=True
			skip_truth = True
			continue
		# skip positions that cannot be typed correctly from no re-typing tool since at least one ALT allele of the sample is missing in population
		if not all(a in alleles + ['0'] for a in true_alleles):
			sys.stderr.write('Skipping position ' + chrom + ':' + position + ' since it is not typable (at least one sample allele is not present in selected samples).\n')
			skip_not_typable += 1
			skip_line=True
			skip_truth=True
			continue
		# skip duplicated or overlapping positions
		if (chrom == prev_chrom) and (int(position) < prev_pos):
			sys.stderr.write('Skipping position ' + chrom + ':' + position + ' since it is duplicated or overlaps a previous variant.\n')
			skip_overlapping += 1
			skip_line=True
			skip_truth = True
			continue
		file.write('\t'.join(vcf_line) + '\n')
	# write truth file if requested
	if skip_line:
		skipped_variants[variant_type] += 1
	if args.truth and not skip_truth:
		truth_file.write('\t'.join(fields[:9] + [fields[sample_id]]) + '\n')
	if prev_chrom != chrom:
		prev_chrom = chrom
		prev_pos = 0
	else:
		prev_pos = max(prev_pos, int(position) + len(fields[3]))
if (args.truth):
	truth_file.close()
for file in panel_files:
	file.close()

sys.stderr.write('\nSkipped non-typable variants:\t' + str(skip_not_typable))
sys.stderr.write('\nSkipped absent variants:\t' + str(skip_absent))
sys.stderr.write('\nSkipped overlapping variants:\t' + str(skip_overlapping) + '\n')
sys.stderr.write('\nSkipped missing alleles:\t' + str(skip_missing) + '\n')

for v in VariantType:
	sys.stderr.write('\nSkipped ' + str(skipped_variants[v]) + ' of ' + str(total_variants[v]) + ' variants in category ' + str(v.name) )

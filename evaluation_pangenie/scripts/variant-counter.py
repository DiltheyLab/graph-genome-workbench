import vcf
import argparse
from variantclassifier import VariantType,  determine_pyvcf_type
from collections import defaultdict

parser = argparse.ArgumentParser(prog='variant-stats.py', description=__doc__)
parser.add_argument('vcf', metavar='VCF', help='VCF to consider.')
parser.add_argument('outfile', metavar='OUT', help='name of the output file (.tsv)')
args = parser.parse_args()

numbers_varianttypes = defaultdict(int)
chromosomes = ["chr" + str(i) for i in range(1,23)] + ['chrX']

for record in vcf.Reader(open(args.vcf, 'r')):
	if not record.CHROM in chromosomes:
		continue
	if record.num_hom_ref == len(record.samples):
		print('Skip position ' + record.CHROM + ':' + str(record.POS) + ' since variant is absent in all samples.')
		continue
	vartype = determine_pyvcf_type(record)
	numbers_varianttypes[vartype] += 1

variants = sorted([v for v in VariantType], key=lambda s: s.name)

# print the table
# header
with open(args.outfile, 'w') as outfile:
	header = ['type', 'total']
	outfile.write('\t'.join(header) + '\n')
	for variant in variants:
		values = [variant.name]
		values.append(str(numbers_varianttypes[variant]))
		outfile.write('\t'.join(values) + '\n')

import vcf
import argparse
from variantclassifier import VariantType, determine_pyvcf_type
from collections import defaultdict

parser = argparse.ArgumentParser(prog='variant-stats.py', description=__doc__)
parser.add_argument('vcf', metavar='VCF', help='VCF to consider.')
parser.add_argument('outdir', metavar='OUT', help='name of the output directory.')
args = parser.parse_args()

numbers_varianttypes = defaultdict(int)
numbers_sample = defaultdict(int)
total_sample = defaultdict(int)

outfiles = {}
first = True

for record in vcf.Reader(filename=args.vcf):
	if first:
		outfiles = {call.sample : open(args.outdir + '/untypable-bubbles/' + call.sample + '-untypable-bubble.bed', 'w') for call in record.samples}
		first = False
	if record.num_hom_ref == len(record.samples):
		print('Skip position ' + record.CHROM + ':' + str(record.POS) + ' since variant is absent in all samples.')
		continue
	vartype = determine_pyvcf_type(record)
	numbers_varianttypes[vartype] += 1

	# store allele_id -> {sample names}
	genotypes = {}
	allele_to_samples = defaultdict(set)
	for call in record.samples:
		genotype_string = call['GT'].replace('/', '|')
		genotype_list = genotype_string.split('|')
		assert len(genotype_list) == 2
		sample_name = call.sample
		allele_to_samples[genotype_list[0]].add(sample_name)
		allele_to_samples[genotype_list[1]].add(sample_name)
		genotypes[sample_name] = genotype_list
	start = int(record.POS) - 1
	end = start + len(record.REF)

	# determine total/untypable variants
	for sample,alleles in genotypes.items():
		# check if variant present in sample
		if alleles[0] not in ['0', '.'] or alleles[1] not in ['0','.']:
			total_sample[(sample, vartype)] += 1

		if '.' in alleles:
			numbers_sample[(sample, vartype)] += 1
			outfiles[sample].write('\t'.join([str(record.CHROM), str(start), str(end)]) + '\n')
			continue

		# check if variant is typable
		if (len(allele_to_samples[alleles[0]]) == 1) and (alleles[0] not in ['0','.']):
			numbers_sample[(sample, vartype)] += 1
			outfiles[sample].write('\t'.join([str(record.CHROM), str(start), str(end)]) + '\n')
			continue
		if (len(allele_to_samples[alleles[1]]) == 1) and (alleles[1] not in ['0', '.']):
			numbers_sample[(sample, vartype)] += 1
			outfiles[sample].write('\t'.join([str(record.CHROM), str(start), str(end)]) + '\n')


#print(numbers_varianttypes)
#print(numbers_sample)
for sample in outfiles:
	outfiles[sample].close()
samples = sorted(set([s[0] for s in total_sample.keys()]))
variants = sorted([v for v in VariantType], key=lambda s: s.name)

# print the table
# header
with open(args.outdir + '/untypable-bubbles.tsv', 'w') as outfile:
	header = ['type'] + [sample + '\t' for sample in samples] + ['total']
	subheader = [''] + ['total','untypable']*len(samples) + ['\t']
	print(' & '.join(['type'] + [sample + ' & ' for sample in samples] + ['total']) + ' \\\ ')
	print(' & '.join(subheader) + ' \\\ ')
	outfile.write('\t'.join(header) + '\n')
	outfile.write('\t'.join(subheader) + '\n')
	for variant in variants:
		values = [variant.name]
		for sample in samples:
			values.append(str(total_sample[(sample, variant)]))
			values.append(str(numbers_sample[(sample, variant)]))
		values.append(str(numbers_varianttypes[variant]))
		print(' & '.join(values) + ' \\\ ')
		outfile.write('\t'.join(values) + '\n')


	


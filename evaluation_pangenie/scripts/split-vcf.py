import sys
import argparse

parser = argparse.ArgumentParser(prog='split-vcf.py', description="cat <vcf-file> | python3 split-vcf.py ")
parser.add_argument('-d', metavar='DIR', default='.', help='output directory.')
parser.add_argument('--include-missing', action='store_true', default=False, help='If set, include positions with missing alleles.')
args = parser.parse_args()

header_lines = []
samples = None
outfiles = None

for line in sys.stdin:
	if line.startswith('##'):
		# store header lines
		header_lines.append(line)
		continue
	if line.startswith('#'):
		# header line
		fields = line.split()
		# determine sample names
		samples = fields[9:]
		outfiles = []
		for s in samples:
			for hap in [0,1]:
				outfile = open(args.d + '/' + s + '-hap' + str(hap) + '.vcf', 'w')
				# write all header/meta information lines
				for h in header_lines:
					outfile.write(h)
				outfile.write('\t'.join(fields[:9] + [s + '_' + str(hap), '\n']))
				outfiles.append(outfile)
		continue
	assert samples is not None
	fields = line.split()
	if 'N' in fields[4] or 'N' in fields[3]:
		# allele sequence contains 'N'
		continue
	untyped = False
	for genotype in fields[9:]:
		# at least one haplotype does not cover region
		if '.' in genotype:
			untyped = True
			break
	if not untyped or args.include_missing:
		# print haplotype of each sample to a separate output file
		variant_columns = fields[:9]
		genotype_columns = fields[9:]
		for i,genotype in enumerate(genotype_columns):
			assert '|' in genotype
			alleles = genotype.split('|')
			for hap in [0,1]:
				outfiles[2*i + hap].write('\t'.join(variant_columns + [alleles[hap], '\n']))

# close files
for file in outfiles:
	file.close()
assert all([f.closed for f in outfiles])

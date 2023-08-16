import pyfaidx
import sys
import argparse

parser = argparse.ArgumentParser(prog='fix-vcf.py', description=__doc__)
parser.add_argument('vcf', metavar='VCF', help='VCF file to be fixed')
parser.add_argument('fasta', metavar='FASTA', help='FASTA reference sequence')
parser.add_argument('--samples', metavar='SAMPLE', help='only keep column for these samples (comma separated)')
args = parser.parse_args()

total_variants = 0
used_variants = 0
keep_columns = []

# read the reference
reference=pyfaidx.Fasta(args.fasta, as_raw=True, sequence_always_upper=True)
for line in open(args.vcf, 'r'):
	if line.startswith('##'):
		# header line
		print(line[:-1])
		continue
	if line.startswith('#'):
		fields = line.split()
		if args.samples is not None:
			keep_columns = [i for i in range(0,9)] + [fields.index(i) for i in args.samples.split(',')]
		else:
			keep_columns = [i for i in range(len(fields))]
		print('\t'.join([fields[i] for i in range(len(fields)) if i in keep_columns]))
		continue
	fields = line.split()
	info = { f.split('=')[0] : f.split('=')[1] for f in fields[7].split(';')}
	chromosome = fields[0]
	position = int(fields[1])-1
	total_variants += 1
	if info['SVTYPE'] == 'DEL' or info['SVTYPE'] == 'INS':
		# get sequence
		alt = info['SEQ']
		if alt == '0':
			sys.stderr.write('Skipping position {pos} since alternative sequence is not given.\n'.format(pos=position))
			continue
		if len(fields) < 10:
			sys.stderr.write('Skipping position {pos} since no genotype is given.\n'.format(pos=position, length=min_len))
			continue
		length = abs(int(info['SVLEN']))
		if length != len(alt):
			sys.stderr.write('Skipping position {pos} since given length contradicts sequence length.\n'.format(pos=position))
			continue
		genotype = fields[-1].split('/')
		if len(genotype) > 1:
			assert len(genotype) == 2
			if genotype[0] == genotype[1]:
				# use phased notation: | instead of /
				fields[-1] = '|'.join(genotype)
			else:
				sys.stderr.write('Skipping position {pos} since it is unphased.\n'.format(pos=position))
				continue
		if info['SVTYPE'] == 'DEL':
			ref_allele = reference[chromosome][position : position+length+1]
			alt_allele = reference[chromosome][position: position + 1]
			if reference[chromosome][position + 1 : position + length + 1] != alt.upper():
				sys.stderr.write('Skipping position {pos} since given sequence does not match reference.\n'.format(pos=position))
				continue
			fields[3] = ref_allele
			fields[4] = alt_allele
		else:
			ref_allele = reference[chromosome][position : position + 1]
			alt_allele = reference[chromosome][position : position +1] + alt.upper()
			assert ref_allele == fields[3]
			fields[3] = ref_allele
			fields[4] = alt_allele
		used_variants += 1
		print('\t'.join([fields[i] for i in range(len(fields)) if i in keep_columns]))
	else:
		sys.stderr.write('Skipping position {pos} since neither INS nor DEL.\n'.format(pos=position))
sys.stderr.write('\n extracted {used} of {total} variants from VCF.'. format(used=used_variants, total=total_variants))


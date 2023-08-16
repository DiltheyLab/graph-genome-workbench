import sys
import argparse

def run_preprocessing(reference_file, vcf_file, length):
	"""
	Adds sequence lengths to header of the VCF file.
	"""
	lengths = {}
	for line in open(reference_file, 'r'):
		fields = line.split()
		lengths[fields[0]] = fields[1]

	for line in open(vcf_file, 'r'):
#		if line.startswith('##contig=<ID'):
#			fields = line.split('=')
#			chrom = fields[2][:-2]
#			print("##contig=<ID=" + chrom + ",length=" + lengths[chrom] + ">")
#			continue

		if not line.startswith('#'):
			fields=line.split()
			ref=fields[3]
			alt=fields[4].split(',')

			alleles = [ref] + alt
			varlen = max([len(a) for a in alleles])
			if varlen >= int(args.l):
				continue

			# skip positions with not containing any IDs of length >= threshold
#			info_fields = {f.split('=')[0] : f.split('=')[1] for f in fields[7].split(';') if '=' in f}
#			assert 'ID' in info_fields
#			allele_ids = info_fields['ID'].split(',')
#			# handle merged IDs
#			all_ids = []
#			for i in allele_ids:
#				for j in i.split(':'):
#					all_ids.append(j)
#			all_ids = list(set(all_ids))
#			varlen = max([int(s.split('-')[-1]) for s in all_ids])
#			if varlen >= int(args.l):
#				continue

		print(line[:-1])

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='gatk-helper.py', description=__doc__)
	subparsers = parser.add_subparsers(dest='subparser_name')
	parser_preprocess = subparsers.add_parser('preprocess', help='prepare vcf to be used as input for paragraph.')
	parser_preprocess.add_argument('fai', metavar='FAI', help='reference fai file.')
	parser_preprocess.add_argument('vcf', metavar='VCF', help='panel VCF file with input variants.')
	parser_preprocess.add_argument('-l', metavar='MAXLENGTH', type=int, default=0, help='max. variant length. Larger variants are excluded.')

	args = parser.parse_args()

	if args.subparser_name == 'preprocess':
		run_preprocessing(args.fai, args.vcf, args.l)


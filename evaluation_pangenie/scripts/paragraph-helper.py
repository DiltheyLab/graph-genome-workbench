import pyfaidx
import sys
import argparse

def run_preprocessing(reference_file, vcf_file):
	"""
	Adds a padding pase (paragraph gives errors instead, even though VCF is valid)
	and modifies the variant position accordingly.
	Positions where a padding base was added are marked as such, such that the
	padding base can be removed again later.
	"""
	reference = pyfaidx.Fasta(reference_file, as_raw=True, sequence_always_upper=True)
	for line in open(vcf_file, 'r'):
		if line.startswith('##'):
			print(line[:-1])
			continue
		if line.startswith('#'):
			fields = line.split()
			print('\t'.join(fields[:8]))
			continue
		else:
			fields=line.split()
			pos=int(fields[1])
			ref=fields[3]
			alt=fields[4].split(',')

			# skip positions with not containing any IDs of length >= threshold
			info_fields = {f.split('=')[0] : f.split('=')[1] for f in fields[7].split(';') if '=' in f}
			assert 'ID' in info_fields
			allele_ids = info_fields['ID'].split(',')
			# handle merged IDs
			all_ids = []
			for i in allele_ids:
				for j in i.split(':'):
					all_ids.append(j)
			all_ids = list(set(all_ids))
			varlen = max([int(s.split('-')[-1]) for s in all_ids])
			if varlen < int(args.l):
				continue

			if not all([len(i) == 1 for i in [ref] + alt]):
				# check if all start with same letter
				if all([a.startswith(ref[0]) for a in alt]):
					print('\t'.join(fields[:8]))
					continue
				else:
					pos=pos-1
					ref=reference[fields[0]][pos-1] + ref
					assert ref == reference[fields[0]][pos-1 : pos-1+len(ref)]
					alt=','.join([reference[fields[0]][pos-1]+a for a in alt])
					fields[1] = str(pos)
					fields[3] = ref
					fields[4] = alt
					fields[2] = 'PADDING'
					print('\t'.join(fields[:8]))
					continue
			else:
				assert len(ref) == 1
				assert all([len(a) == 1 for a in alt])
				print('\t'.join(fields[:8]))

def run_postprocessing(vcf_file):
	"""
	Removes the padding base in regions marked with
	PADDING in ID field of the input VCF.
	Respective positions are increased by 1.
	"""
	for line in open(vcf_file, 'r'):
		if line.startswith('#'):
			print(line[:-1])
			continue
		else:
			fields=line.split()
			base_added = ('PADDING' in fields[2])
			if not base_added:
				print(line[:-1])
				continue
			# there is a padding base, remove it
			pos=int(fields[1])
			ref=fields[3]
			alts=fields[4].split(',');
			# remove base
			fields[1] = str(pos + 1)
			fields[2] = '.'
			fields[3]=ref[1:]
			alts = ','.join([a[1:] for a in alts])
			fields[4] = alts
			print('\t'.join(fields))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='paragraph-helper.py', description=__doc__)
	subparsers = parser.add_subparsers(dest='subparser_name')
	parser_preprocess = subparsers.add_parser('preprocess', help='prepare vcf to be used as input for paragraph.')
	parser_preprocess.add_argument('reference', metavar='REFERENCE', help='reference sequence in fasta-format.')
	parser_preprocess.add_argument('vcf', metavar='VCF', help='panel VCF file with input variants.')
	parser_preprocess.add_argument('-l', metavar='MINLENGTH', default=0, help='only consider variants of at least this length.')

	parser_postprocess = subparsers.add_parser('postprocess', help='fix the paragraph output vcf.')
	parser_postprocess.add_argument('vcf', metavar='VCF', help='paragraph VCF output.')
	args = parser.parse_args()

	if args.subparser_name == 'preprocess':
		run_preprocessing(args.reference, args.vcf)
	if args.subparser_name == 'postprocess':
		run_postprocessing(args.vcf)


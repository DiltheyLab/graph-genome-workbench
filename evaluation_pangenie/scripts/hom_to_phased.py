#!/usr/bin/env python3

## this script changes notation of hom genotypes from 1/1 -> 1|1
## also checks if chromosomes are phased in a single block

import sys


def main():
	current_chrom = None
	current_block = None
	for line in sys.stdin:
		if line.startswith('#'):
			print(line[:-1])
			continue
		fields = line.split()
		# assume single-sample vcf
		assert len(fields) == 10
		if current_chrom != fields[0]:
			current_chrom = fields[0]
			current_block = None
		sample = fields[9].split(':')
		genotype = sample[0]
		if (len(sample) > 1) and sample[-1] != '.':
			block_id = sample[1]
			if current_block is None:
				current_block = block_id
			if (block_id != current_block):
				sys.stderr.write('skip position ' + fields[0] + ':' + fields[1] + ' since it is located in a different block.\n')
				continue
		alleles = [int(i) for i in genotype.replace('/', '|').split('|')]
		if len(set(alleles)) == 1:
			sample[0] = '|'.join([str(i) for i in alleles])
		fields[9] = ':'.join(sample)
		print('\t'.join(fields))


if __name__ == "__main__":
	main()

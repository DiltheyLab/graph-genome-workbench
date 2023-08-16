#!/usr/bin/env python3

## this script removes all variants from the given vcf for which the distance
## to previous and next variant is at least the given threshold

import sys
import argparse

def main():
	# read provided vcf to be filtered
	parser = argparse.ArgumentParser(prog='remove_clusters.py', description=__doc__)
	parser.add_argument('threshold', metavar='THRESHOLD', help='remove all variants closer than THRESHOLD.')
	args = parser.parse_args()

	current_cluster = []
	prev_pos = 0
	prev_chrom = None

	for line in sys.stdin:
		if line[0] == '#':
			# header line, skip
			print(line[:-1])
			continue
		
		columns = line.split('\t')
		position = int(columns[1])
		chrom = columns[0]

		if ((prev_chrom == chrom) and ((position - prev_pos) < int(args.threshold))) or (prev_chrom == None):
			current_cluster.append(line)
		else:
			# variant does not belong to current cluster
			# if current cluster contains only a single variant, print it
			if len(current_cluster) == 1:
				print(current_cluster[0][:-1])
			current_cluster = [line]
		if prev_chrom == chrom:
			prev_pos = max(prev_pos, position + len(columns[3]))
		else:
			prev_pos = position + len(columns[3])
		prev_chrom = chrom
	
	# check if to print last cluster
	if len(current_cluster) == 1:
		print(current_cluster[0][:-1])


if __name__ == "__main__":
	main()

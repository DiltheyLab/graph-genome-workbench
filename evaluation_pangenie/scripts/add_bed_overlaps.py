import sys, argparse
from collections import defaultdict
import pandas as pd

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='add_bed_overlaps.py', description=__doc__)
	parser.add_argument('-t', '--table', required=True, help='Table to which annotations shall be added.')
	parser.add_argument('-i', '--intersections', required=True, nargs='+', help='bedtools intersect output for two VCF files (keeping both entries)')
	parser.add_argument('-o', '--output', required=True, help='output table with added columns.')
	parser.add_argument('-n', '--name', default='intersection', nargs='+', help='name of the callset that was intersected.')
	args = parser.parse_args()

	df = pd.read_csv(args.table, sep='\t')
	# add table indicating overlap
	for filename, name in zip(args.intersections, args.name):
		callset_name = name
		intersections = defaultdict(list)
		for line in open(filename, 'r'):
			fields = line.split()
			id_first = fields[2]
			id_second = fields[8] + ':' + fields[9] + '-' + fields[10]
			intersections[id_first].append(id_second)
		in_callset = []
		id_callset = []
		for var in df['ID']:
			if var in intersections:
				in_callset.append(True)
				id_callset.append(','.join(intersections[var]))
			else:
				in_callset.append(False)
				id_callset.append('nan')
		df['in_' + callset_name] = in_callset
		df['ID_' + callset_name] = id_callset
	df.to_csv(args.output, sep='\t', index=False, na_rep='nan')

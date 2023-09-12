## The difference betwwen annotate.py and annotate-ids.py is that the former anotates a VCF file based on IDs of other file and the latter takes the ID from the info fields and set it to the ID within the same VCF file
## Additionally, annotate-ids.py is designed for biallelic VCF files

import sys

# read VCF and store IDs of all alleles
for line in sys.stdin:
	if line.startswith('#'):
		print(line[:-1])
		continue
	fields = line.split()
	alleles = fields[4].split(',')
	assert len(alleles) == 1 # only biallelic files
	info_field = {f.split('=')[0] : f.split('=')[1] for f in fields[7].split(';') if '=' in f}
	assert 'ID' in info_field
	ids = info_field['ID'].split(',')
	assert len(ids) == len(alleles)
	fields[2] = ids[0] # since len(ids) == 1
	print('\t'.join(fields))		


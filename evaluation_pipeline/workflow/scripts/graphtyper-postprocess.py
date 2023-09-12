import sys
import gzip

input_vcf = sys.argv[1]

id_to_alleles = {}

for line in gzip.open(input_vcf, 'rt'):
	if line.startswith('#'):
		continue
	fields = line.split()
	info_fields = {s.split('=')[0] : s.split('=')[1] for s in fields[7].split(';') if '=' in s}
	assert 'ID' in info_fields
	alleles = fields[4].split(',')
	# assume biallelic VCF
	assert len(alleles) == 1
	id_to_alleles[info_fields['ID']] = [fields[3], alleles[0]]


for line in sys.stdin:
	if line.startswith('#'):
		print(line.strip())
		continue
	fields = line.strip().split()
	info_fields = {s.split('=')[0] : s.split('=')[1] for s in fields[7].split(';') if '=' in s}
	assert 'SVMODEL' in info_fields
	if info_fields['SVMODEL'] != 'AGGREGATED':
		continue
	assert 'OLD_VARIANT_ID' in info_fields
	var_id = info_fields['OLD_VARIANT_ID']
	assert var_id in id_to_alleles
	fields[3] = id_to_alleles[var_id][0]
	fields[4] = id_to_alleles[var_id][1]
	print('\t'.join(fields))
	

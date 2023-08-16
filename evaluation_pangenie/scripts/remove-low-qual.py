import sys

threshold = float(sys.argv[1])

for line in sys.stdin:
	if line.startswith('#'):
		print(line[:-1])
		continue
	fields = line.split()
	assert len(fields) > 9
	format = fields[8].split(':')
	assert 'GQ' in format
	assert 'GT' in format
	index_gq = format.index('GQ')
	index_gt = format.index('GT')
	for i in range(len(fields[9:])):
		genotyping = fields[i+9].split(':')
		quality = genotyping[index_gq]
		genotype = genotyping[index_gt]
		if quality != '.':
			if float(quality) < threshold:
				genotyping[index_gt] = './.'
		fields[i+9] = ':'.join(genotyping)
	print('\t'.join(fields))

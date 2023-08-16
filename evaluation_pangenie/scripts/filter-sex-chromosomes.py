import sys
from collections import defaultdict

samples = []

male_samples = []
female_samples = []
set_missing = defaultdict(lambda: 0)
skip_no_alt = 0

total = 0
written = 0

for line in open(sys.argv[1], 'r'):
	fields = line.strip().split()
	if fields[1] == 'm':
		male_samples.append(fields[0])
	else:
		female_samples.append(fields[0])


for line in sys.stdin:
	fields = line.strip().split()
	if line.startswith('##'):
		print(line[:-1])
		continue
	if line.startswith('#'):
		samples = fields[9:]
		print(line[:-1])
		continue
	total += 1
	nonref = False
	if 'X' in fields[0]:
		# check if males have 1|1 or 0|0 genotypes
		for i in range(len(samples)):
			sample = samples[i]
			genotype = fields[9+i]
			if genotype in ['0|1', '1|0', '1/0', '0/1'] and sample in male_samples:
				sys.stderr.write('Wrong genotype for sample ' + sample + ' at position ' + fields[0] + ':' + fields[1] + '.\n')
				fields[9+i] = '.|.'
				set_missing[sample] += 1
			if '1' in fields[9+i]:
				nonref = True
	elif 'Y' in fields[0]:
		# check if female samples have 0|0 genotypes
		for i in range(len(samples)):
			sample = samples[i]
			genotype = fields[9+i]
			if not genotype in ['0|0', '0/0'] and sample in female_samples:
				sys.stderr.write('Wrong genotype for sample ' + sample + ' at position ' + fields[0] + ':' + fields[1] + '.\n')
				fields[9+i] = '.|.'
				set_missing[sample] += 1
			if '1' in fields[9+i]:
				nonref = True
	else:
		nonref=True

	if nonref:
		print('\t'.join(fields))
		written += 1
	else:
		sys.stderr.write('No alt alleles at position ' + fields[0] + ':' + fields[1] + '.\n')
		skip_no_alt += 1

for sample in samples:
	sys.stderr.write('Set ' + str(set_missing[sample]) + ' of sample ' + sample + ' to missing.\n')

sys.stderr.write('Skip no alt positions ' + str(skip_no_alt) + '.\n')

sys.stderr.write('kept ' + str(written) + ' of ' + str(total) + ' variants.\n')

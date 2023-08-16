import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+')
args = parser.parse_args()

total_runtime = 0
for filename in args.files:
	with open(filename, 'r') as file:
		for line in file:
			if 'User time (seconds)' in line or 'System time (seconds)' in line:
				fields = line.split()
				total_runtime += float(fields[-1])
print('total: ' + str(total_runtime))

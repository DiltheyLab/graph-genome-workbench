import sys
from variantclassifier import VariantType, determine_variant_from_line

import argparse

parser = argparse.ArgumentParser(prog='extract-varianttype.py', description=__doc__)
parser.add_argument('vartype', metavar='TYPE', help='variant type. Possible options: snp|small-insertion|small-deletion|midsize-insertion|midsize-deletion|large-insertion|large-deletion|indel|sv|small|midsize|large|all')
args = parser.parse_args()

string_to_type = {
    "indel": [VariantType.snp, VariantType.small_insertion, VariantType.small_deletion, VariantType.small_complex, VariantType.midsize_insertion, VariantType.midsize_deletion, VariantType.midsize_complex],
	"sv": [VariantType.large_insertion, VariantType.large_deletion, VariantType.large_complex],
    "small": [VariantType.snp, VariantType.small_insertion, VariantType.small_deletion, VariantType.small_complex],
	"midsize" : [VariantType.midsize_insertion, VariantType.midsize_deletion, VariantType.midsize_complex],
	"large": [VariantType.large_insertion, VariantType.large_deletion, VariantType.large_complex]
}

assert args.vartype in string_to_type
type_to_extract = string_to_type[args.vartype]

total = 0
written = 0

for line in sys.stdin:
	if line.startswith('#'):
		print(line.strip())
		continue
	vartype = determine_variant_from_line(line)
	fields = line.strip().split()
	if vartype in type_to_extract:
		print(line.strip())
		written += 1
	total += 1
sys.stderr.write('Extracted ' + str(written) + ' of ' + str(total) + ' variants.\n')
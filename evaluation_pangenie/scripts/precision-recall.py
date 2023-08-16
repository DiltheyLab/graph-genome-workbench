import sys
import argparse
import vcf
from enum import Enum
from multiprocessing import Pool

class VariantType(Enum):
	snp = 0
	small_insertion = 1
	small_deletion = 2
	small_complex = 3
	midsize_insertion = 4
	midsize_deletion = 5
	midsize_complex = 6
	large_insertion = 7
	large_deletion = 8
	large_complex = 9

def determine_variant_type(record):
	"""
	Determines the variant type.
	"""
	alleles = [record.REF] + record.ALT
	varlen = max([len(a) for a in alleles])

	if record.is_snp:
		return VariantType.snp

	is_deletion = record.var_subtype == 'del'
	is_insertion = record.var_subtype == 'ins'

	if varlen < 20:
		if is_insertion:
			return VariantType.small_insertion
		if is_deletion:
			return VariantType.small_deletion
		return VariantType.small_complex

	if varlen >= 20 and varlen <= 50:
		if is_insertion:
			return VariantType.midsize_insertion
		if is_deletion:
			return VariantType.midsize_deletion
		return VariantType.midsize_complex

	if varlen > 50:
		if is_insertion:
			return VariantType.large_insertion
		if is_deletion:
			return VariantType.large_deletion
		return VariantType.large_complex


def read_vcf(filename):
	"""
	reads the vcf and returns a list with counts
	of each variant type.
	"""
	counts = [0 for vartype in VariantType]
	for record in vcf.Reader(open(filename, 'rb')):
		vartype = determine_variant_type(record)
		counts[vartype.value] += 1
	return counts


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='precision-recall.py', description=__doc__)
	parser.add_argument('tp_baseline', metavar='TPBASELINE', help='VCF with true positive baseline variants.')
	parser.add_argument('tp_callset', metavar='TPCALLSET', help='VCF with true positive callset variants.')
	parser.add_argument('fp', metavar='FP', help='VCF with false positive variants.')
	parser.add_argument('fn', metavar='FN', help='VCF with false negative variants.')
	args = parser.parse_args()

	pool = Pool(4)

	# counts of tp, fp and fn per variant type
	tp_baseline_res = pool.apply_async(read_vcf, [args.tp_baseline])
	tp_callset_res = pool.apply_async(read_vcf, [args.tp_callset])
	fp_res = pool.apply_async(read_vcf, [args.fp])
	fn_res = pool.apply_async(read_vcf, [args.fn])

	tp_baseline_variants = tp_baseline_res.get()
	tp_callset_variants = tp_callset_res.get()
	fp_variants = fp_res.get()
	fn_variants = fn_res.get()

	print('\t'.join(['variant_type', 'TP_baseline', 'TP_callset', 'FP', 'FN', 'precision', 'recall', 'F-score']))

	# compute precision and recall for each variant type
	for vartype in VariantType:
		tp_baseline = tp_baseline_variants[vartype.value]
		tp_callset = tp_callset_variants[vartype.value]
		fp = fp_variants[vartype.value]
		fn = fn_variants[vartype.value]

		recall = tp_baseline / max(tp_baseline + fn, 1)
		precision = tp_callset / max(tp_callset + fp, 1)
		f_score = ( (2*precision*recall) / (precision + recall) ) if (precision + recall) > 0.0 else 0

		print('\t'.join([vartype.name, str(tp_baseline), str(tp_callset), str(fp), str(fn), str(precision), str(recall), str(f_score)]))

	# results for all variant types combined
	tp_baseline = sum(tp_baseline_variants)
	tp_callset = sum(tp_callset_variants)
	fp = sum(fp_variants)
	fn = sum(fn_variants)

	recall = tp_baseline / max(tp_baseline + fn, 1)
	precision = tp_callset / max(tp_callset + fp, 1)
	f_score = ( (2*precision*recall) / (precision + recall) ) if (precision + recall) > 0.0 else 0

	print('\t'.join(['all', str(tp_baseline), str(tp_callset), str(fp), str(fn), str(precision), str(recall), str(f_score)]))


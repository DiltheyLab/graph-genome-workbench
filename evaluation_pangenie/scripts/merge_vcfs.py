#!/usr/bin/python

import sys
import argparse
import re
import pyfaidx
from collections import defaultdict
from variantclassifier import VariantType, determine_type_from_ids

class Genotype:
	def __init__(self, alleles, is_phased):
		self._alleles = alleles
		self._phased = is_phased

	def __str__(self):
		allele_str = []
		for allele in self._alleles:
			if allele == -1:
				allele_str.append('.')
			else:
				allele_str.append(allele)
		if self._phased:
			return '|'.join([str(i) for i in allele_str])
		else:
			return '/'.join([str(i) for i in allele_str])

	def get_alleles(self):
		return self._alleles

	def get_ploidy(self):
		return len(self._alleles)

	def is_phased(self):
		return self._phased

	def __eq__(self, other):
		if self._phased != other._phased:
			return False
		if self._phased:
			return self._alleles == other._alleles
		else:
			return sorted(self._alleles) == sorted(other._alleles)

	def is_hom_ref(self):
		return all([int(a) == 0 for a in self._alleles])


class Variant:

	def __init__(self, samples, start, ref_allele, alt_alleles, genotypes, ploidy, var_ids):
		assert len(samples) == len(genotypes)
		assert len(alt_alleles) == len(var_ids)
		self._samples = samples
		self._start = start
		self._end = start + len(ref_allele)
		self._ref_allele = ref_allele
		self._alt_alleles = alt_alleles
		self._genotypes = genotypes
		self._ploidy = ploidy
		self._id = var_ids

	def get_samples(self):
		return self._samples

	def get_start(self):
		return self._start

	def get_end(self):
		return self._end

	def get_ref(self):
		return self._ref_allele

	def get_alt(self):
		return self._alt_alleles

	def get_ploidy(self):
		return self._ploidy

	def get_id(self, allele_index):
		if allele_index == 0:
			assert(False)
		else:
			return self._id[allele_index-1]

	def get_ids(self):
		return self._id

	def get_allele(self, index):
		if index == 0:
			return self._ref_allele
		else:
			return self._alt_alleles[index - 1]

	def get_genotype_of(self, sample):
		if sample in self._samples:
			return self._genotypes[self._samples.index(sample)]
		else:
			return Genotype([0] * self._ploidy, True)

	def is_absent(self):
		alleles = []
		for g in self._genotypes:
			for a in g.get_alleles():
				alleles.append(a)
		return all(x == 0 for x in alleles)

	def __eq__(self, other):
		if self._start != other._start:
			return False
		if self._end != other._end:
			return False
		if self._samples != other._samples:
			return False
		if self._ref_allele != other._ref_allele:
			return False
		if self._alt_alleles != other._alt_alleles:
			return False
		if self._genotypes != other._genotypes:
			return False
		if self._ploidy != other._ploidy:
			return False
		if self._id != other._id:
			return False
		return True

	def __lt__(self, other):
		return self._start < other._start

	def __le__(self, other):
		return self._start <= other._start

	def __gt__(self,other):
		return self._start > other._start

	def __ge__(self, other):
		return self._start >= other._start

	def __str__(self):
		sample_str = '[' + ','.join([str(s) for s in self._samples]) + ']'
		start_str = str(self._start)
		end_str = str(self._end)
		alt_str = ','.join(self._alt_alleles)
		genotype_str = '[' + ','.join([str(g) for g in self._genotypes]) + ']'
		id_str = '[' + ','.join(self._id) + ']'
		return 'Variant(' + ','.join([sample_str, start_str, end_str, self._ref_allele, alt_str, genotype_str, id_str]) + ')'


def insert_variant(sequence, variant, allele_index, offset):
	start = variant.get_start() - offset
	end = variant.get_end() - offset
	included = sequence[0:start] + variant.get_allele(allele_index) + sequence[end:]
	offset -= len(included) - len(sequence)
#	print(sequence, start, variant.get_allele(allele_index), included)
	return included, offset

def combine_variants(variants, reference, chromosome):
	if all(v.is_absent() for v in variants):
		return None
	# determine start, end
	start = min([v.get_start() for v in variants if not v.is_absent()])
	end = max([v.get_end() for v in variants if not v.is_absent()])
#	var_id = ';'.join(set([v.get_id() for v in variants if not v.is_absent()]))

	# determine samples
	samples = set([])
	for v in variants:
		for s in v.get_samples():
			samples.add(s)

	samples = sorted(list(samples))

	# determine reference sequence
	ref_allele = reference[chromosome][start-1 : end-1]
	
	# ploidy
	ploidy = variants[0].get_ploidy()
	
	genotypes = []
	new_alleles = defaultdict(lambda:-1)
	new_alleles[ref_allele] = 0
	new_ids = {}
	allele_index = 1

	for sample in samples:
		genotype_alleles = []
		for p in range(0,ploidy):
			haplotype = ref_allele
			ids = ''
			offset = start
			prev_end = 0
			undefined_allele = False
			for v in variants:
				allele = v.get_genotype_of(sample).get_alleles()[p]
				# if undefined allele, the combined allele is undefined as well
				if allele == -1:
					undefined_allele = True
					continue
				if allele != 0:
					if prev_end > v.get_start():
						sys.stderr.write('Two overlapping variants on the same haplotype at: ' + chromosome + ':' + str(v.get_start()) + ' ' + ','.join(v.get_ids()) + ' . Set allele to missing.\n')
						undefined_allele = True
						continue
					haplotype, offset = insert_variant(haplotype, v, allele, offset)
					if ids != '':
						ids += ':'
					ids += v.get_id(allele)
					prev_end = v.get_end()
			if undefined_allele:
				genotype_alleles.append(-1)
				
				continue
			if new_alleles[haplotype] == -1:
				new_alleles[haplotype] = allele_index
				new_ids[haplotype] = ids
				allele_index += 1
			index = new_alleles[haplotype]
			genotype_alleles.append(index)
		genotypes.append(Genotype(genotype_alleles, True))
				
	alt_alleles = [''] * len(new_alleles)
	var_ids = [''] * len(new_alleles)
	for allele, index in new_alleles.items():
		alt_alleles[index] = allele
		if index != 0:
			# reference allele does not have an ID
			var_ids[index] = new_ids[allele]

	return Variant(list(samples), start, ref_allele, alt_alleles[1:], genotypes, ploidy, var_ids[1:])	
			

def genotype_from_string(gt_string):
	is_phased = False
	alleles = []
	if '|' in gt_string:
		is_phased = True
		alleles = []
		for allele in gt_string.split('|'):
				if allele != '.':
					alleles.append(int(allele))
				else:
					alleles.append(-1)
	elif '/' in gt_string:
		is_phased = False
		for allele in gt_string.split('/'):
				if allele != '.':
					alleles.append(int(allele))
				else:
					alleles.append(-1)
	else:
#		assert len(gt_string) == 1
		is_phased = True
		alleles = [int(gt_string)] if gt_string != '.' else [-1]
	return Genotype(alleles, is_phased)


def parse_line(samples, line, id_in_info=False):
	variants = []
	fields = line.split()
	chrom = fields[0]
	start = int(fields[1])
	ref_allele = fields[3]
	alt_alleles = fields[4].split(',')
	ids = fields[2].split(';')
	if id_in_info:
		# pav vcfs have ids in INFO column...
		info = { i.split('=')[0] : i.split('=')[1] for i in fields[7].split(';') if "=" in i}
		assert 'ID' in info
		ids = info['ID'].split(',')
	for i,s in enumerate(samples):
		# get samples genotype
		genotype = genotype_from_string(fields[9 + i])
		if not genotype.is_phased():
			raise Exception('Line: ' + line + ' contains an unphased genotype.')
		variants.append(Variant([s], start, ref_allele, alt_alleles, [genotype], genotype.get_ploidy(), ids))
	return chrom, variants

def combine_sorted_lists(list1, list2):
	combined = []
	size1 = len(list1)
	size2 = len(list2)
	i,j = 0,0
	while i < size1 and j < size2:
		if list1[i] < list2[j]:
			combined.append(list1[i])
			i += 1
		else:
			combined.append(list2[j])
			j += 1
	return combined + list1[i:] + list2[j:]
	

def print_header(samples):
	print('##fileformat=VCFv4.2')
	print('##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs per ALT allele.\">')
	print('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">')
	print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples))

# return: skipped_conflict, total_written
def print_cluster(variants, samples, chrom, reference):
	# skipped due to conflict
	skipped_conflict = defaultdict(int)
	# variants written to output
	total_written = defaultdict(int)

	if len(variants) == 0:
		return skipped_conflict, total_written
	combined = combine_variants(variants, reference, chrom)
	if combined is not None:
		start = str(combined.get_start())
		ref_allele = combined.get_ref()
		alt_alleles = ','.join(combined.get_alt())
		variant_ids = ','.join(combined.get_ids())
		if len(alt_alleles) < 1:
			sys.stderr.write('Skipping variant at position ' + chrom + ':' + start + ' '  +  variant_ids + ' since no alternative alleles are present after merging.\n')
			return skipped_conflict, total_written
		if any(c not in 'CAGTcagt,' for c in ref_allele) or any(c not in 'CAGTcagt,' for c in alt_alleles):
			sys.stderr.write('Skipping variant at position ' + chrom + ':' + start + ' '  +  variant_ids + ' since invalid character in allele sequence found.\n')
			return skipped_conflict, total_written
		genotypes = '\t'.join([str(combined.get_genotype_of(s)) for s in samples])
		outline = '\t'.join([chrom, start, '.', ref_allele, alt_alleles, '.', 'PASS', 'ID=' + variant_ids, 'GT', genotypes])
		print(outline)
		# print how many variants are left
		allele_ids = set([])
		for id_string in combined.get_ids():
				for id in id_string.split(':'):
					allele_ids.add(id)
		types = determine_type_from_ids(allele_ids)
		for v_type in types:
			total_written[v_type] += 1
		return skipped_conflict, total_written
		
	else:
		# count how many variants are in cluster
		allele_ids = set([])
		for v in variants:
			for id_string in v.get_ids():
				for id in id_string.split(':'):
					allele_ids.add(id)
		types = determine_type_from_ids(allele_ids)
		for v_type in types:
			skipped_conflict[v_type] += 1	
		return skipped_conflict, total_written



def run_merge(reference, vcfs, ploidy, output_chromosomes):
	# statistics per variant type
	total_input = defaultdict(int)
	skipped_conflict = defaultdict(int)
	total_written = defaultdict(int)

	ref = pyfaidx.Fasta(reference, as_raw=True, sequence_always_upper=True)
	# store chromosome --> [list of variants]
	variants_per_chromosome = defaultdict(list)
	all_samples = []
	# read vcf-files
	for filename in vcfs:
		vcf = open(filename, 'r')
		samples = None
		for line in vcf:
			if line.startswith('##'):
				continue # change here if we want to include header lines
			if line.startswith('#'):
				# determine sample names
				fields = line.split()
				if samples is not None:
					raise Exception('File ' + filename + ' is not a valid VCF file.')
				if len(fields) < 10:
					raise Exception('File ' + filename + ' does not contain any samples.')
				samples = fields[9:]
				continue
			chromosome, variants = parse_line(samples, line, True)
			# include new variants if present in at least one sample (ignore positions with only 0|0 genotypes)
			if not all(v.is_absent() for v in variants):
				variants_per_chromosome[chromosome].extend(variants)
#			variants_per_chromosome[chromosome] = combine_sorted_lists(variants_per_chromosome[chromosome], variants)
		all_samples += samples
		vcf.close()

	# sort variants
	for chrom in variants_per_chromosome:
		variants_per_chromosome[chrom] = sorted(variants_per_chromosome[chrom])

	# write VCF header
	print_header(all_samples)
	# write variants
	chromosomes = sorted(variants_per_chromosome.keys())
	if len(output_chromosomes) > 0:
		chromosomes = sorted(output_chromosomes)
	for chromosome in chromosomes:
		variants = variants_per_chromosome[chromosome]
		current_cluster = []
		prev_end = 0
		all_ids = set([])
		for v in variants:
			# count alleles
			for id_string in v.get_ids():
				for id in id_string.split(':'):
					all_ids.add(id)
			if v.get_start() >= prev_end:
				skipped, written = print_cluster(current_cluster, all_samples, chromosome, ref)
				for k,val in skipped.items():
					skipped_conflict[k] += val
				for k,val in written.items():
					total_written[k] += val
				current_cluster = []
			current_cluster.append(v)
			prev_end = max(v.get_end(), prev_end)
		# print last cluster
		skipped, written = print_cluster(current_cluster, all_samples, chromosome, ref)
		for k,val in skipped.items():
			skipped_conflict[k] += val
		for k,val in written.items():
			total_written[k] += val
		types = determine_type_from_ids(all_ids)
		for v_type in types:
			total_input[v_type] += 1

	# print statistics
	for vartype in VariantType:
		sys.stderr.write(vartype.name + ': Skipped due to conflicts: ' + str(skipped_conflict[vartype]) + ', Written: ' + str(total_written[vartype]) + ', Total input: ' + str(total_input[vartype]) + '\n')
			

def combine_haplotypes(variants):
	alleles = []
	for v in variants:
		samples = v.get_samples()
		assert len(samples) == 1
		genotype = v.get_genotype_of(samples[0])
		alleles += genotype.get_alleles()
	return Genotype(alleles, True)
		

def run_combine(vcf, sample):
	with open(vcf, 'r') as outfile:
		nr_samples = None
		for line in outfile:
			if line.startswith('##'):
				print(line[:-1])
				continue
			if line.startswith('#'):
				samples = line.split('\t')[9:]
				print('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + [sample]))	
				continue
			chromosome, variants = parse_line(samples,line,True)
			assert len(variants) == len(samples)
			fields = line.split('\t')
			assert len(fields) > 9
			fields[8] = 'GT'
			fields[9] = str(combine_haplotypes(variants))
			print('\t'.join(fields[:10]))

def run_combine_columns(vcf, samples):
	sample_to_haplotype = defaultdict(list)
	for line in open(samples, 'r'):
		fields = line.split()
		sample_to_haplotype[fields[0]] = fields[1:]
	column_to_index = {}
	with open(vcf, 'r') as outfile:
		nr_samples = None
		for line in outfile:
			if line.startswith('##'):
				print(line[:-1])
				continue
			if line.startswith('#'):
				columns = line[:-1].split('\t')[9:]
				for index,c in enumerate(columns):
					column_to_index[c] = index
				samples = [s for s in sample_to_haplotype.keys()]
				print('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + samples))	
				continue
			vcf_samples = []
			for s,v in sample_to_haplotype.items():
				vcf_samples.extend(v)
			chromosome, variants = parse_line(vcf_samples,line,True)
			fields = line.split('\t')[:9]
			fields[8] = 'GT'
			for sample,haps in sample_to_haplotype.items():
				fields.append(str(combine_haplotypes([variants[column_to_index[h]] for h in haps])))
			print('\t'.join(fields))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='merge_vcfs.py', description=__doc__)
	subparsers = parser.add_subparsers(dest="subparser_name")
	parser_merge = subparsers.add_parser('merge', help='merge single-sample VCF files into one multi-sample VCF file.')
	parser_merge.add_argument('-r', metavar='reference', required=True, help='reference sequence (FASTA-format).')
	parser_merge.add_argument('-vcf', metavar='VCF', required=True, nargs='+', help='VCF-file(s) with variants to merge.')
	parser_merge.add_argument('-ploidy', metavar='PLOIDY', required=True, type=int, help='ploidy of the samples.')
	parser_merge.add_argument('-chromosomes', metavar='CHROMOSOMES', default='', help='comma separated list of chromosomes. Only output variants on these chromosomes.')
	parser_combine = subparsers.add_parser('combine', help='combine single haplotype columns into a genotype column.')
	parser_combine.add_argument('-vcf', metavar='VCF', required=True, help='VCF-file with one column per haplotype.' )
	parser_combine.add_argument('-sample', metavar='SAMPLE', help='name of the sample to appear in final VCF-file.', default='sample')
	parser_combine = subparsers.add_parser('combine_columns', help='combine single haplotype columns into diploid genotype columns.')
	parser_combine.add_argument('-vcf', metavar='VCF', required=True, help='VCF-file with one column per haplotype.' )
	parser_combine.add_argument('-samples', metavar='SAMPLES', help='file containing one line per sample to be merged (sample_name,h0,h1)', required=True)
	args = parser.parse_args()

	if args.subparser_name == 'merge':
		chromosomes = []
		if args.chromosomes != '':
			chromosomes = args.chromosomes.split(',')
		run_merge(args.r, args.vcf, args.ploidy, chromosomes)
	if args.subparser_name == 'combine':
		run_combine(args.vcf, args.sample)
	if args.subparser_name == 'combine_columns':
		run_combine_columns(args.vcf, args.samples)

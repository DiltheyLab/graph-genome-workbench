#!/usr/bin/python

import argparse
import vcf
import sys
from collections import namedtuple
from collections import defaultdict
from variantclassifier import determine_variant_length, determine_variant_class

AlleleStats = namedtuple('AlleleStats','af ac an')
GenotypeStats = namedtuple('GenotypeStats', 'heterozygosity het total')

def compute_allele_statistics(record):
	"""
	Compute allele related statistics.
	"""
	an = 0
	ac = 0
	for genotype in record.samples:
		alleles = genotype['GT'].replace('/', '|').split('|')
		assert 1 <= len(alleles) <= 2
		if len(alleles) == 1:
			# haploid genotype
			alleles.append(alleles[0])
		for a in alleles:
			if a == '.':
				continue
			assert int(a) in [0,1]
			an += 1
			ac += int(a)
	if an < 1:
		assert ac < 1
	af = ac / max(1.0, float(an))
	return AlleleStats(str(af), str(ac), str(an))


def read_uk(record):
	return str(record.INFO['UK'])
	

def compute_genotype_statistics(record, qualities=None):
	"""
	Compute genotype related statistics.
	"""
	counts = defaultdict(int)
	het_genotypes = 0
	total_genotypes = 0
	gqs = [record.samples[i]['GQ'] for i in range(len(record.samples))] if qualities is not None else [None]*len(record.samples)
	for genotype, quality in zip(record.samples, gqs):
		alleles = genotype['GT'].replace('/', '|').split('|')
		assert 1 <= len(alleles) <= 2
		if len(alleles) == 1:
			# haploid genotype
			alleles.append(alleles[0])
		if not '.' in alleles:
			alleles = [int(a) for a in alleles]
			total_genotypes += 1
			if sum(alleles) == 1:
				assert 0 in alleles
				assert 1 in alleles
				het_genotypes += 1
			# read GQ
			if qualities is not None:
				for q in qualities:
					if int(quality) >= q:
						counts[q] += 1
	genotype_stats = GenotypeStats( str(het_genotypes / max(1.0, float(total_genotypes))), str(het_genotypes), str(total_genotypes))
	return genotype_stats, counts


parser = argparse.ArgumentParser(prog='collect-vcf-stats.py', description=__doc__)
parser.add_argument('panel', metavar='panel', help='biallelic panel variants.')
parser.add_argument('pangenie', metavar='pangenie', help='PanGenie biallelic genotyped variants.')
args = parser.parse_args()

panel_stats = {}

for variant in vcf.Reader(open(args.panel,'r')):
	# require bi-allelic vcf with IDs
	assert len(variant.ALT) == 1
	var_id = variant.INFO['ID']
	assert len(var_id) == 1
	allele_stats = compute_allele_statistics(variant)
	panel_stats[var_id[0]] = allele_stats

pangenie_stats = {}
quals = [0,200]

for variant in vcf.Reader(open(args.pangenie,'r')):
	assert len(variant.ALT) == 1
	var_id = variant.INFO['ID']
	assert len(var_id) == 1
	allele_stats = compute_allele_statistics(variant)
	genotype_stats, counts = compute_genotype_statistics(variant, quals)
	uk = read_uk(variant)
	pangenie_stats[var_id[0]] = [allele_stats, genotype_stats, uk, counts]


# print stats for all IDs in genotypes VCF
header = [ 	'variant_id',
		'panel_allele_freq',
		'panel_alternative_alleles',
		'panel_total_alleles',
		'pangenie_allele_freq',
		'pangenie_alternative_alleles',
		'pangenie_total_alleles',
		'pangenie_heterozygosity',
		'pangenie_heterozygous_genotypes',
		'pangenie_total_genotypes',
		'pangenie_unique_kmers'
	]

for q in quals:
	header.append('pangenie_GQ>=' + str(q))

print('\t'.join(header))

for var_id in pangenie_stats:
	if not var_id in panel_stats:
		continue
		
	line = [	var_id,
			panel_stats[var_id].af,
			panel_stats[var_id].ac,
			panel_stats[var_id].an,

			pangenie_stats[var_id][0].af,
			pangenie_stats[var_id][0].ac,
			pangenie_stats[var_id][0].an,
			pangenie_stats[var_id][1].heterozygosity,
			pangenie_stats[var_id][1].het,
			pangenie_stats[var_id][1].total,
			pangenie_stats[var_id][2]
		]
		
	# add counts for GQs
	for q in quals:
		line.append(str(pangenie_stats[var_id][3][q]))
	assert len(line) == len(header)
	print('\t'.join(line))

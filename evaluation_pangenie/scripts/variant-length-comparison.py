import matplotlib
matplotlib.use('pdf')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
import vcf

from variantclassifier import determine_variant_length, determine_variant_class, determine_class_from_line

#plt.style.use('seaborn-deep')


def extract_variant_from_line(line):
	vartype = determine_class_from_line(line)
		
	splitted = line.split()
	alleles = [splitted[3]] + [s for s in splitted[4].split(',')]
	varlen = max([len(a) for a in alleles])
	if (vartype == 'INS') or (vartype == 'DEL'):
		# insertions and deletions always have one base more
		varlen -= 1
	assert varlen >= 1
	return vartype, varlen
	

def extract_AF(line):
	fields = line.split()
	info_fields = { i.split('=')[0] : i.split('=')[1].strip() for i in fields[7].split(';') if '=' in i}
	assert "AF" in info_fields
	return float(info_fields['AF'])


def extract_length(record, is_symbolic):
	fields = record.split()
	if is_symbolic:
		info = fields[7]
		assert "SVLEN" in info
		info = {i.split('=')[0] : i.split('=')[1] for i in info.split(';') if "SVLEN" in i} 
		svlengths = [int(i) for i in info['SVLEN'].split(',')]
		return svlengths
	else:
		alleles = [fields[3]] + fields[4].split(',')
		svlengths = [len(i) for i in alleles]
		return svlengths

parser = argparse.ArgumentParser(prog='variant-length-comparison.py')
parser.add_argument('gnomad', metavar="GNOMAD", help="VCF containing Gnomad variants")
parser.add_argument('assemblies', metavar="ASSEMBLIES", help="VCF containing assembly variants.")
parser.add_argument('-outprefix', metavar="OUTPREFIX", help="prefix of the output files.", default="output")
args = parser.parse_args()

af_threshold = 0.05

gnomad_hists = {'all': [], 'ins':[], 'del': []}
assembly_hists = {'all': [], 'ins':[], 'del': []}

# read gnomad vcf
for line in open(args.gnomad, 'r'):
	if line[0] == '#':
		continue
	if line[0] == 'Y':
		continue
	f = line.split('\t')
	varlen = max(extract_length(line, True))
	gnomad_hists['all'].append(varlen)

	# extract allele frequency from AF
	allele_freq = extract_AF(line)
	if allele_freq < af_threshold:
		continue

	if f[4] == "<DEL>":
		gnomad_hists['del'].append(-varlen)
	elif f[4] == "<INS>":
		gnomad_hists['ins'].append(varlen)
	elif f[4] == "<DUP>":
		gnomad_hists['ins'].append(varlen)
	elif "<CN" in f[4]:
		gnomad_hists['ins'].append(varlen)

total_count = 0
# write assembly vcf with variants >= 50 only
with open(args.outprefix + "/assembly-svs.vcf", "w") as outvcf:
	for line in open(args.assemblies, 'r'):
		if line[0] == '#':
			outvcf.write(line)
			continue
		variant_type, varlen = extract_variant_from_line(line)
		if variant_type in ['INS', 'DEL'] and varlen >= 50:
			outvcf.write(line)
			assembly_hists['all'].append(varlen)
			total_count += 1

complex_count = 0
del_count = 0
ins_count = 0
skipped = 0
for record in vcf.Reader(open(args.outprefix + "/assembly-svs.vcf", "r")):
	variant_type = determine_variant_class(record)
	varlen = determine_variant_length(record)
	assert varlen >= 50

	# determine allele frequency
	assert len(record.INFO['AF']) == 1
	allele_freq = float(record.INFO['AF'][0])
	assert 0.0 <= allele_freq <= 1.0

	if allele_freq < af_threshold:
		skipped += 1
		continue

	if variant_type == 'DEL':
		assembly_hists['del'].append(-varlen)
		del_count += 1
	elif variant_type == "INS":
		assembly_hists['ins'].append(varlen)
		ins_count += 1
	else:
		complex_count += 1

assert total_count == ins_count + del_count + complex_count + skipped
print('deletions', del_count, 'insertions', ins_count, 'all', total_count)

with PdfPages(args.outprefix + "/length-histogram.pdf") as pdf:
	for svtype in ['all', 'ins', 'del']:
		all_numbers = [assembly_hists[svtype], gnomad_hists[svtype]]
		bins = np.logspace(1.6,8,200) if svtype != 'del' else -1*np.logspace(1.6,8,200)[::-1]
		for numbers, label, color in zip(all_numbers, ['PanGenie', 'gnomAD'], ['red','blue']):
			y,binEdges=np.histogram(numbers,bins=bins)
			bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
			plt.plot(bincenters,y,'-', label=label, color=color)
		plt.title(svtype)
		plt.xscale('symlog')
		plt.yscale('symlog')
		plt.xlabel('variant length')
		plt.ylabel('count')
		plt.legend(loc='upper right')
		pdf.savefig()
		plt.close()		

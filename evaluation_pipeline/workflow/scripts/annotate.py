import sys
import argparse
from collections import defaultdict
import gzip

parser = argparse.ArgumentParser(prog='annotate.py', description=__doc__)
parser.add_argument('vcf', metavar='ANNOTATIONS', help='VCF to take annotations from.')
args = parser.parse_args()

# store (chrom, pos) -> allele -> id
pos_to_alleles = defaultdict(dict)

alleles_read = 0
alleles_assigned = 0
unknown_alleles = 0
alleles_skipped = 0

# read VCF and store IDs of all alleles
#for line in open(args.vcf, 'rt'):
for line in gzip.open(args.vcf, 'rt'):
    if line.startswith('#'):
        if line.startswith('##INFO=<ID=ID'):
            id_line = line
        continue
    fields = line.split()
    chrom = fields[0]
    pos = int(fields[1])
    alleles = fields[4].split(',')
    info_field = {f.split('=')[0] : f.split('=')[1] for f in fields[7].split(';') if '=' in f}
    assert 'ID' in info_field
    ids = info_field['ID'].split(',')
    assert len(ids) == len(alleles)
    for allele, id in zip(alleles, ids):
        assert not allele in pos_to_alleles[(chrom, pos, fields[3])]
        pos_to_alleles[(chrom, pos, fields[3])][allele] = id
        alleles_read += 1
sys.stderr.write('Read ' + str(alleles_read) + ' alleles/ids from ' + args.vcf + '.\n')


# read other VCF and add IDs to matching alleles
index = 0
for line in sys.stdin:
    if line.startswith('##'):
        there_is_not_id_line = False if line.startswith('##INFO=<ID=ID') else True            
        print(line[:-1]) 
        continue
    if line.startswith('#'):
        if there_is_not_id_line:
            print(id_line[:-1]) ## inserting ID INFO tag in the header
        print(line[:-1]) 
        continue
    fields = line.strip().split()
    chrom = fields[0]
    pos = int(fields[1])
    alleles = fields[4].split(',')
    info_field = {f.split('=')[0] : f.split('=')[1] for f in fields[7].split(';') if '=' in f}
    ids = []
    if '*' in alleles:
        # skip line
        alleles_skipped += len(alleles)
        continue
    for i,allele in enumerate(alleles):
        if allele in pos_to_alleles[(chrom, pos, fields[3])]:
            id = pos_to_alleles[(chrom, pos, fields[3])][allele]
            ids.append(id)
            alleles_assigned += 1
        else:
            ids.append('-'.join(['allele', chrom, str(pos), fields[3], allele]))
            unknown_alleles += 1
            #sys.stderr.write(str(pos) + '\n')
            #sys.stderr.write(str(ids) + '\n\n')
            #continue
    info_field['ID'] = ','.join(ids)
    if len(ids) == 1:
        fields[2] = ids[0]
    fields[7] = ';'.join([f + '=' + info_field[f] for f in info_field])
    index += 1
    print('\t'.join(fields))
sys.stderr.write('Assigned ' + str(alleles_assigned) + ' allele ids to variants in input VCF.\n')
sys.stderr.write(str(unknown_alleles) + ' alleles were unknown.\n')
sys.stderr.write(str(alleles_skipped) + ' alleles were skipped due to *.\n')

import sys
from collections import defaultdict

total_lines = 0
skipped_lines = 0
for line in sys.stdin:
    if line.startswith('#'):
        print(line[:-1])
        continue
    fields = line.split()
    ref = fields[3]
    alts = fields[4]
    total_lines += 1
    
    # Check if the REF allele contains letters other than "ATGC"
    if any(letter not in 'ATGC' for letter in ref):
        skipped_lines += 1
        continue  # Skip this record
    
    # Check if any ALT allele contains letters other than "ATGC"
    if any(any(letter not in 'ATGC' for letter in alt) for alt in alts):
        skipped_lines += 1
        continue  # Skip this record
    
    print(line[:-1])

sys.stderr.write('Filtered ' + str(skipped_lines) + ' of ' + str(total_lines) + ' from input vcf because no ATGC.\n')	

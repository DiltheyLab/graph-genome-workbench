import sys

for line in sys.stdin:
    if line.startswith('##'):
        print(line.strip())
        continue
    print('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency.">')
    print(line.strip())

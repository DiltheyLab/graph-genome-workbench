import sys

lines_to_append = ['##contig=<ID=chr1,length=248956422>',
                    '##contig=<ID=chr2,length=242193529>',
                    '##contig=<ID=chr3,length=198295559>',
                    '##contig=<ID=chr4,length=190214555>',
                    '##contig=<ID=chr5,length=181538259>',
                    '##contig=<ID=chr6,length=170805979>',
                    '##contig=<ID=chr7,length=159345973>',
                    '##contig=<ID=chr8,length=145138636>',
                    '##contig=<ID=chr9,length=138394717>',
                    '##contig=<ID=chr10,length=133797422>',
                    '##contig=<ID=chr11,length=135086622>',
                    '##contig=<ID=chr12,length=133275309>',
                    '##contig=<ID=chr13,length=114364328>',
                    '##contig=<ID=chr14,length=107043718>',
                    '##contig=<ID=chr15,length=101991189>',
                    '##contig=<ID=chr16,length=90338345>',
                    '##contig=<ID=chr17,length=83257441>',
                    '##contig=<ID=chr18,length=80373285>',
                    '##contig=<ID=chr19,length=58617616>',
                    '##contig=<ID=chr20,length=64444167>',
                    '##contig=<ID=chr21,length=46709983>',
                    '##contig=<ID=chr22,length=50818468>',
                    '##contig=<ID=chrX,length=156040895>']

def check_lines_in_header(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        not any(line.strip() in header_lines for line in lines)

def append_lines_to_vcf_header(file):
    
    with open(file, 'r+') as f:
        lines = f.readlines()
        if not any(line.strip() in header_lines for line in lines):
            f.seek(0, 0)  # Move the file pointer to the beginning
            for line in reversed(header_lines):
                f.write(line + '\n')
            f.writelines(lines)
            print(f"Lines appended to {file}")

# Usage example
vcf_files = ['file1.vcf', 'file2.vcf', 'file3.vcf']
append_lines_to_vcf_header(vcf_files)

parser = argparse.ArgumentParser(prog='append_lines_header.py', description=__doc__)
parser.add_argument('vcf', metavar='input VCF', help='VCF to append header lines to, if missing')
args = parser.parse_args()



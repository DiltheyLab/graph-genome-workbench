#!/bin/sh

fasta=$1
vcf=$2
output=$3

python3 ../scripts/paragraph-helper.py preprocess $fasta $vcf -l 20 | bgzip > $output

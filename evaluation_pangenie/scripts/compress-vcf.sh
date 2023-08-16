#!/bin/sh

vcffile=$1
output=$2

bgzip -c $vcffile > $output
tabix -p vcf $output

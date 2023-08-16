#!/bin/sh

platypus=$1
bam=$2
fasta=$3
output=$4
vcf=$5
log=$6

export C_INCLUDE_PATH=/home/ebler/miniconda3/envs/platypus-install/include
export LIBRARY_PATH=/home/ebler/miniconda3/envs/platypus-install/lib
export LD_LIBRARY_PATH=/home/ebler/miniconda3/envs/platypus-install/lib

(/usr/bin/time -v python2 $platypus callVariants --bamFiles=$bam --refFile=$fasta --output=$output  --source=$vcf --minPosterior=0 --getVariantsFromBAMs=0) &> $log

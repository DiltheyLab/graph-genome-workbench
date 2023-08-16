#!/bin/sh

input_bam=$1
output_bam=$2
output_metrics=$3
log=$4

/usr/bin/time -v picard MarkDuplicates INPUT=$input_bam OUTPUT=$output_bam METRICS_FILE=$output_metrics &> $log
samtools index $output_bam

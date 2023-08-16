#!/bin/sh

region=$1
computed=$2
bed=$3
output=$4

if [ "$region" == "genotyping-non-repetitive" ]
then
	bedtools intersect -f 1.0 -header -a $computed -b $bed | uniq > $output
else
	bedtools intersect -header -a $computed -b $bed | uniq > $output
fi

#!/bin/sh

coverage=$1
fraction=$2
reads1=$3
reads2=$4
out1=$5
out2=$6

cov=$(<$coverage)

p=$(awk -v frac=$fraction -v c=$cov 'BEGIN { printf "%.2f", frac/c }')

printf "p=%f\n" $p

[ "$fraction" -gt $cov ] && echo "Not enough coverage" && exit 1

seqtk sample $reads1 $p > $out1
seqtk sample $reads2 $p > $out2

exit 0

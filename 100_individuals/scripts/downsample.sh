#!/bin/bash

#pipeline is: split vcf into the proper chromosomes
#Need to accont for the unknown data

# dir=$1
for ind in {1..100}; do
  n=`wc -l ${ind}.2L.dat | awk '{print $1}'`
  n=$((n*5/100))
  echo $ind
  `cat <(cat <(shuf ${ind}.2L.dat | head -n $n | awk '{print $1,$2}') <(shuf ${ind}.2L.dat | head -n $n | awk '{print $1,$3}') | head -n $n | sort -nuk 1) <(awk '{print $1,"."}' ${ind}.2L.dat ) | sort -nuk 1 > ${ind}.2L.downsampled.dat`
  sed -i '1s/.*/Pos obs/' ${ind}.2L.downsampled.dat

done

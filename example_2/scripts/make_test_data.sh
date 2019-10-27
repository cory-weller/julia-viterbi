#!/bin/bash

examples=("1" "2" "3")
chromosomes=("2L" "2R" "3L" "3R" "X")

for example in ${examples[@]}; do
	for chromosome in ${chromosomes[@]}; do
		echo "Doing $example for $chromosome"
		head -n 2 ${example}_observations_${chromosome}.dat > ${example}_observations_${chromosome}_summary.dat
		head -n 2 ${example}_${chromosome}_true_founders.dat > ${example}_${chromosome}_true_founders_summary.dat
		awk 'NR % 10 == 0' ${example}_observations_${chromosome}.dat >> ${example}_observations_${chromosome}_summary.dat
		awk 'NR % 10 == 0' ${example}_${chromosome}_true_founders.dat >> ${example}_${chromosome}_true_founders_summary.dat	
	done
done

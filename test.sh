#!/bin/bash
# Badri Adhikari, 10/15/2015

for i in `seq 1 23`; do
	echo "Running job for Chromosome ${i} at 1MB.."
	./chromosome3D.pl -if "./input/chr${i}_1mb_matrix.txt" -o "../output/chr${i}_1mb" &> "../output/chr${i}_1mb.log" &
done

for i in `seq 1 23`; do
	echo "Running job for Chromosome ${i} at 500KB.."
	./chromosome3D.pl -if "./input/chr${i}_500kb_matrix.txt" -o "../output/chr${i}_500kb" &> "../output/chr${i}_500kb.log" &
done
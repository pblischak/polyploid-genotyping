#!/bin/bash

for f in *.fq.gz
do
	PREFIX=$(echo $f | awk -F'.' '{print $1}')
	bwa mem -t 2 ../../Betula_concat_reference.fasta $f | \
	samtools view -bSu - | \
	samtools sort -O bam -o ../bam/$PREFIX.sorted.bam
done

#!/bin/bash

samtools mpileup -I -f ../../Betula_concat_reference.fasta \
  -l ../../filtered30-variants.txt $(ls *.sortedRG.bam) \
  -o ../../filtered30-pubescens.pileup

#!/bin/bash

# Script for filtering a VCF file and converting to read count matrices.
# The script requires VCFtools and Python.

# To see the options for the Python script, type: ./read-counts-from-vcf.py -h

# This is the command that was used to filter the data for the Andropogon gerardii data analysis.
vcftools --gzvcf McAllister.Miller.all.mergedRefGuidedSNPs.vcf.gz \
         --remove-indv Blank:C28DEACXX:6:250193414 \
         --max-alleles 2 --min-alleles 2 --thin 10000 \
         --minDP 5 --max-missing 0.5 --remove-indels \
         --remove-filtered-all --recode \
         --stdout | ./read-counts-from-vcf.py --ad 2 --minDP 5

#!/bin/bash

# Script for filtering a VCF file and converting to read count matrices.
# The script requires VCFtools and Perl.

# The Perl script is pretty simple and is run with the following options
# passed to it (in order):

# perl read-counts-from-vcf.pl <total-reads-file> <alt-reads-file> <allele-depth-position> <min-depth-filter>

#   total-reads-file       -- the name to be given to the total read count file.
#   alt-reads-file         -- the name to be given to the alternative read count file.
#   allele-depth-position  -- the position of the allele depth information in the genotype field for the VCF file.
#   min-depth-filter       -- the minimum number of reads for including read count data.

# This is the command that was used to filter the data for the Andropogon gerardii data analysis.
vcftools --gzvcf McAllister.Miller.all.mergedRefGuidedSNPs.vcf.gz \
         --remove-indv Blank:C28DEACXX:6:250193414 \
         --max-alleles 2 --min-alleles 2 --thin 10000 \
         --minDP 5 --max-missing 0.5 --remove-indels \
         --remove-filtered-all --recode \
         --stdout | perl read-counts-from-vcf.pl test-tot.txt test-ref.txt 2 5

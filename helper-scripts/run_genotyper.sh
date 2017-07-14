#!/bin/bash

# Use the GATK UnifiedGenotyper to make a VCF file
# for all BAM files.

java -jar ~/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
	-T UnifiedGenotyper \
	-R ../../Betula_concat_reference.fasta \
	$(for f in *sortedRG.bam; do printf "%sI $f " '-'; done) \
	-o ../pendula-ug.vcf

# Test on a single BAM file
#java -jar ~/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
#	-T HaplotypeCaller \
#	-R ../../Betula_concat_reference.fasta \
#	-I 1147x_CTCTCTAG.sortedRG.bam \
#	-o test.vcf

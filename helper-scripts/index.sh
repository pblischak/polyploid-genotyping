#!/bin/bash

# Create sequence dictionary for use with GATK
java -jar ~/picard-tools-2.2.1/picard.jar CreateSequenceDictionary R=Betula_concat_reference.fasta O=Betula_concat_reference.dict

# Create SAMtools index for use with GATK
samtools faidx Betula_concat_reference.fasta

# Create BWA index for read mapping
bwa index Betula_concat_reference.fasta

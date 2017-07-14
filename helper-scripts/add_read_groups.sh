#!/bin/bash

for f in *.bam
do
	PREFIX=$(echo $f | awk -F'.' '{print $1}')
	java -jar ~/picard-tools-2.2.1/picard.jar AddOrReplaceReadGroups \
		I=$f \
		O=$PREFIX.sortedRG.bam \
		SORT_ORDER=coordinate \
		RGID=pendula \
		RGLB=pendula \
		RGPL=illumina \
		RGSM=$PREFIX \
		RGPU=pendula \
		CREATE_INDEX=True
done
#!/bin/bash


FILENAME=BNvs392.new.merged.vqsr.pass.se.
GATK=/home/nkong/software/GenomeAnalysisTK-3.7-0-gcfedb67/GenomeAnalysisTK.jar
REF=/data/epfl/fellay/aborghes/databases/old_fasta/human_g1k_v37_decoy.fa
SnpSift=/mnt/data1/HIV_SC/nkong/snpEff/SnpSift.jar


java -Xmx8g -jar $GATK \
	-R $REF \
	-T VariantFiltration \
	--variant $FILENAME.vcf \
	-o $FILENAME.tmp1.vcf \
	--filterExpreslssion "AC > 0" \
	--filterName "ALTER"



cat $FILENAME.vcf | java -jar $SnpSift filter " F(ILTER = 'ALTER')" > $FILENAME.filt.vcf


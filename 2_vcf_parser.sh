#!/bin/bash

##Script uses vcf_parser to split lines where multiple variants with different AFs are present at the same locus  

#module add $PATH/vcf_parser;

#FILENAME is initialised in all_script
#if you are not using all_scripts then decomment the next line and write the name of the file

FILENAME=BNvs392.new.merged.vqsr.pass.se.
DIR=/mnt/data1/HIV_SC/nkong/


#PARSING - old way (it takes a long time; BCFtools is much faster) 
vcf_parser --split $FILENAME.vcf > $FILENAME.split.vcf


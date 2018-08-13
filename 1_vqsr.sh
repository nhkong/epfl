#!/bin/bash


# Script filters a vcf according to per sample (not per batch) DP and GQ at a given site; 

SnpSift=/home/nkong/snpEff/SnpSift.jar

DIR=/mnt/data1/HIV_SC/nkong/

FILENAME=BNvs392.new.merged.vqsr.pass.se


cat $DIR/$FILENAME.vcf | java -jar $SnpSift filter "( GEN[*].DP > 10 ) & ( GEN[*].GQ > 20 )" > $DIR/$FILENAME.vcf`

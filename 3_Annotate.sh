#!/bin/bash

##Script performs multiple functional annotation commands: 

# 1 Add Variant Type (SNP/MNP/INS/DEL/MIXED) in the INFO field with SnpSift
# 2 Annotate using dbNSFP 
# 3 Annotate Add LOF field 
# 4 Rename files 
# 5 Remove intermediate files 

#FILENAME is initialised in all_script
#if you are not using all_scripts then decomment the next line and write the name of the file
FILENAME=BNvs392.new.merged.vqsr.pass.se.

SnpSift=/home/nkong/snpEff/SnpSift.jar
SNPEFF=/home/nkong/snpEff/snpEff/snpEff.jar

extension="split.vcf"

# 1 

java -Xmx10g -jar $SnpSift varType $FILENAME.$extension > $FILENAME.$extension.vt
extension="$extension".vt
# 2

#added Interpro_domain:domain or conserved site on which the Qvariant locates. 

#added LRT_pred:LRT prediction, D(eleterious), N(eutral) or U(nknown)

#added: 1000Gp1_EUR_AF: in order to also compare only to the european descendent sample

#additional info: http://varianttools.sourceforge.net/Annotation/DbNSFP 
#if you want to add something you should just change the next line and the script number 4 (header and fields)
java -Xmx10g -jar $SnpSift dbnsfp -v -db /data/epfl/fellay/databases/dbNSFP/2.9.3/dbNSFP2.9.3.txt.gz -f rs_dbSNP147,refcodon,aaref,aaalt,aapos,Ancestral_allele,ExAC_AF,ExAC_AC,clinvar_trait,clinvar_rs,CADD_raw,CADD_raw_rankscore,CADD_phred,MetaSVM_rankscore,MetaLR_rankscore,SIFT_converted_rankscore,SIFT_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_pred,LRT_converted_rankscore,MutationTaster_converted_rankscore,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,MutPred_Top5features,SiPhy_29way_logOdds_rankscore,GERP++_RS_rankscore,phyloP46way_primate_rankscore,phastCons46way_primate_rankscore,1000Gp1_AF,1000Gp1_EUR_AF,ESP6500_EA_AF,ESP6500_AA_AF,ARIC5606_EA_AF,ARIC5606_AA_AF $FILENAME.$extension > $FILENAME.$extension.dbNSFP
extension="$extension".dbNSFP

# 3 

java -Xmx10g -jar $SNPEFF ann -lof GRCh37.75 -c /data/epfl/fellay/aborghes/software/snpEff/snpEff/snpEff.config $FILENAME.$extension > $FILENAME.$extension.LOF
extension="$extension".LOF

# 4

file=$FILENAME.$extension
mv -i "${file}" "${file/vcf.vt.dbNSFP.LOF/vt.dbNSFP.LOF.vcf}"

# 5

rm *.vt
rm *.dbNSFP

# END`

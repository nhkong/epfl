```#!/bin/bash

#FILENAME and sampleIDs are initialised in all_script
#if you are not using all_scripts then decomment the 2 next lines and write the name of the file and samples
FILENAME=BNvs392.new.merged.vqsr.pass.se
sampleIDs="PIBB_2501   PIBB_2502       PIBB_2503"
extension=split.vt.dbNSFP.LOF
OUTPUT=$FILENAME
SNPSIFT=/home/nkong/snpEff/snpEff/SnpSift.jar


#1- extract fields with snpsift
#from a VCF file to a TXT, tab separated file (that can be easily loaded to XLS and R)
#-e to specify how are empty fields
#VARTYPE: {SNP/INS/DEL...}   dbNSFP_ESP6500_EA_AF: alternative allele frequency in european/american samples
#dbNSFP_1000Gp1_AF: filter variants that are probably damaging     dbNSFP_SIFT_pred: predicts if damaging or telerated
#GEN[*].GT: genotype subfields from all samples      LOF[*].GENE: estimate if loss of function
#added all fields annotated from script 3
#IF YOU WANT TO ADD ANNOTAION YOU SHOULD JUST 1)CHANGE SCRIPT 3 AND 2)PUT THE NAME IN FIELDS AND 3)IN THE HEADER AT THE END OF THIS SCRIPT
#"GEN[*].GT" "LOF[*].GENE" ANN should be the 3 last ones (and GENE_NAME the last one before other annotation)

fields=(CHROM POS ID REF ALT AC AF VARTYPE dbNSFP_rs_dbSNP147 dbNSFP_refcodon dbNSFP_aaref dbNSFP_aaalt dbNSFP_aapos dbNSFP_Ancestral_allele dbNSFP_clinvar_trait dbNSFP_clinvar_rs dbNSFP_M_CAP_score dbNSFP_M_CAP_rankscore dbNSFP_REVEL_score dbNSFP_REVEL_rankscore dbNSFP_MutPred_score dbNSFP_MutPred_rankscore dbNSFP_SIFT_converted_rankscore dbNSFP_Polyphen2_HVAR_score dbNSFP_MutationTaster_converted_rankscore dbNSFP_LRT_converted_rankscore dbNSFP_MetaLR_rankscore dbNSFP_MetaSVM_rankscore dbNSFP_phastCons46way_primate_rankscore dbNSFP_phyloP46way_primate_rankscore dbNSFP_GERP_RS_rankscore dbNSFP_SiPhy_29way_logOdds_rankscore dbNSFP_ESP6500_EA_AF dbNSFP_ESP6500_AA_AF dbNSFP_1000Gp1_AF dbNSFP_1000Gp1_EUR_AF dbNSFP_CADD_raw dbNSFP_CADD_raw_rankscore dbNSFP_CADD_phred "dbNSFP_Polyphen2_HVAR_pred" dbNSFP_LRT_pred dbNSFP_SIFT_pred dbNSFP_M_CAP_pred dbNSFP_MutPred_Top5features "GEN[*].GT" "GEN[*].DP" dbNSFP_ExAC_AF "NMD[*].GENE" "LOF[*].GENE" ANN)

#declares header
#needs to be in the same order as FIELDS
echo "CHR:POS:REF:ALT GENE ID REF ALT AC AF VARTYPE rs_dbSNP147 refcodon aaref aaalt aapos Ancestral_allele clinvar_trait clinvar_rs M-CAP_score M-CAP_rankscore REVEL_score REVEL_rankscore MutPred_score MutPred_rankscore SIFT_converted_rankscore Polyphen2_score MutationTaster_converted_rankscore LRT_converted_rankscore LR_rankscore RadialSVM_rankscore phastCons_primate_rankscore phyloP_primate_rankscore GERP_rankscore SiPhy_rankscore EA_AF AA_AF 1kG_AF 1kG_EUR_AF CADD_raw CADD_raw_rankscore CADD_phred Polyphen2 LRT SIFT M-CAP MutPred "$sampleIDs" "$sampleIDs" ExACv1_AF NMD LOF CLASS IMPACT GENE_NAME HGVS.c HGVS.p Hom_exac_v2 Hom_gnomad ExACv2 GNOMAD INHOUSE mcap_v1" > 15.tmp
cat 15.tmp | tr " " "\t" > 15.1.tmp

# where duplicate column names are present (e.g. originating from sample names "$sampleIDs" "$sampleIDs" called twice), Pandas.py script adds a suffix to make column name unique 
# e.g. NI_101 NI_101 become NI_101 NI_101.1
# output is file 16.tmp
python Pandas.py

#Changes the GERP++_RS to GERP because extracting fields with SNIPSIFT can have some problems if none alphanumeric characters (i.e: ++)
cat $FILENAME.$extension.vcf | sed "s/dbNSFP_GERP__/dbNSFP_GERP/g" > $FILENAME.$extension.GERP.vcf

#removes variants were alt allele is *
#awk -F '\t' '{if ($5 ~ /\*/) print$0}' $FILENAME.$extension.GERP.vcf > $FILENAME.$extension.GERP.nodel.vcf

java -jar $SNPSIFT extractFields $FILENAME.$extension.GERP.vcf ${fields[@]} > $FILENAME.$extension.tmp


#2- getting rid of header line
#+ inverts the tail (everything but the 2 first lines)
#needs to put in tmp because > opens before tail
tail -n +2 $FILENAME.$extension.tmp > 0.tmp


#3- removing tab between chr number and variant position
cut -f1,2 0.tmp | sed 's/\t/:/g' > 1.tmp


#4- geti some of the columns
lastF=$(awk -F '\t' '{print NF}' 0.tmp | head -n 1)
cut -f3-$(($lastF-1)) 0.tmp > 2.tmp

#garde seulement annotation, impact et gene name de la fin= keeps only annotation, impact and gene name of the end
#integration: keep also HGVD.p and HGVD.c notations for the variants 
cut -f"$lastF" 0.tmp | sed 's/|/\t/g' | cut -f2-11 > 3.tmp
cut -f1 3.tmp | sed 's/^$/NotAv/g' > 3.1.tmp
cut -f2 3.tmp | sed 's/^$/NotAv/g' > 3.2.tmp
cut -f3 3.tmp | sed 's/^$/NF/g' > 3.3.tmp
cut -f10 3.tmp | sed 's/^$/NotAv/g' > 3.4.tmp
cut -f11 3.tmp | sed 's/^$/NotAv/g' > 3.5.tmp
paste 3.1.tmp 3.2.tmp 3.3.tmp 3.4.tmp 3.5.tmp > 4.tmp

paste 1.tmp 2.tmp 4.tmp > 5.tmp


#5- change "Chr pos rsID REF ALT etc." in "Chr:pos:REF:ALT rsID etc." 
awk -F "\t" 'BEGIN{OFS=":"}{print $1,$3,$4}' 5.tmp > 6.tmp

#keeps everything from the second field
cut -f2- 5.tmp > 7.tmp

#generates the file that will be annotated later with variants in format Chr:pos:REF:ALT
#adds a GENE column after Chr:pos:REF:ALT
paste -d "\t" 6.tmp 3.3.tmp 7.tmp > 8.tmp


#6- Annotate with AF in ExAc, in the SHCS, and in our in-house exome database
#Please note: a double annotation is performed for each database (2 files per each database): one with Chr:pos:REF:ALT and one with Chr:pos:ALT:REF
#This is because some variants are annotated in the inverted way 'Chr:pos:ALT:REF'. The two columns are then merged to one with a py script (see later)   
#Please also note that this part of the script requires the Annotate_maf.sh script (Samira's script) in the current directory and the databases   
#all the database files are stored in the databases folder in /data/epfl/fellay/aborghes/databases/AFs and are currently copied to the current folder for use in the script and then removed to avoid duplication 
folderAF="/data/epfl/fellay/aborghes/databases/AFs/"
./Annotate_maf.sh 8.tmp "$folderAF"gnomad.exomes.r2.0.1.sites.PASS.split.Hom.txt > 9.tmp
./Annotate_maf.sh 9.tmp "$folderAF"gnomad.genomes.r2.0.1.sites.allchrom.PASS.split.Hom.txt > 10.tmp
./Annotate_maf.sh 10.tmp "$folderAF"gnomad.exomes.r2.0.1.sites.PASS.split.AF.1.txt > 11.tmp
#./Annotate_maf.sh 11.tmp "$folderAF"gnomad.exomes.r2.0.1.sites.PASS.split.AF.2.txt > 12.tmp
./Annotate_maf.sh 11.tmp "$folderAF"gnomad.genomes.r2.0.1.sites.Allchrom.PASS.split.AF.1.txt > 13.tmp
#./Annotate_maf.sh 13.tmp "$folderAF"gnomad.genomes.r2.0.1.sites.Allchrom.PASS.split.AF.2.txt > 14.tmp
./Annotate_maf.sh 13.tmp "$folderAF"HIV.merged.vqsr.pass.split.AF.1.txt > 25.tmp
#./Annotate_maf.sh 25.tmp "$folderAF"HIV.merged.vqsr.pass.split.AF.2.txt > 26.tmp
#./Annotate_doit.sh 25.tmp 31.tmp > 32.tmp
./Annotate_maf.sh 25.tmp "$folderAF"mcap_v1.txt > 26.tmp 


#7- Add header
#needs to be in the same order as FIELDS
cat 16.tmp 26.tmp > 27.tmp

# Because from the last version of GATK some variants that are GT ./. have no other info (e.g. the DP is not reported)
# the script coverage.py doesn't work because automatically SnpSift extractfields substitutes the absebce of the field with NotAv
# Therefore, I have to change NotAv in 0, and I do it in a manual way: 
# Above, I chage the output filename from $OUTPUT.snpSIFTextr_anno.txt in 27.tmp
# And here: 
lastF=$(awk -F '\t' '{print NF}' 27.tmp | head -n 1)
cut -f1-220 27.tmp > 28.tmp
cut -f221-396 27.tmp | sed 's/NotAv/0/g' > 29.tmp
cut -f397-"$lastF" 27.tmp > 30.tmp
paste -d "\t" 28.tmp 29.tmp 30.tmp > $OUTPUT.snpSIFTextr_anno.txt````

#!/bin/bash

# awk to filter out Ref variants  
# awk to filter out Synonymous variants 
# awk according to AFs 

#FILENAME is initialised in all_script
#if you are not using all_scripts then DECOMMENT the next line and write the name of the file
FILENAME="BNvs392.new.merged.vqsr.pass.se"

# filter out variants were called genotypes have DP (read depth) < 10 (kept if at least one genotype called as ALT has DP > 10) 
# Uses py script written by John Wilson on March 2017 
./coverageFilter.py < $FILENAME.snpSIFTextr_anno.txt > $FILENAME.HighCov.txt


# only awk variants that are alt allele in at least one sample 
awk -F '\t' '{if ($0 ~ /0\/1/ || $0 ~ /1\/0/ || $0 ~ /1\/1/ || NR == 1) print$0}' $FILENAME.HighCov.txt > $FILENAME.alt.txt

# awk to filter out synonymous variants
awk -F '\t' '{if ($0 !~ /synonymous_variant/ || NR == 1) print$0}' $FILENAME.alt.txt  > $FILENAME.alt.NonSyn.txt


# awk to filter out according to allele frequency in ExAc and in the inhouse databases
#head is an array containing all the header in order to find the column corresponding to a header (ex: LOF)

declare -A head
eval "$(head -n 1 $FILENAME.alt.NonSyn.txt | awk -F '\t' '{for(i=1;i<=NF;i++) printf "head[%s]=%d ", $i, i }')"


awk -F '\t' ' $"'${head[ExACv2]}'" < 0.01 || $"'${head[ExACv2]}'" ~ /NF/ || NR == 1 ' $FILENAME.alt.NonSyn.txt > $FILENAME.alt.NonSyn.exac.txt


awk -F '\t' ' $"'${head[GNOMAD]}'" < 0.05 || $"'${head[GNOMAD]}'" ~ /NF/ || NR == 1 ' $FILENAME.alt.NonSyn.exac.txt > $FILENAME.alt.NonSyn.exac.gnomad.txt


awk -F '\t' ' $"'${head[INHOUSE]}'" < 0.05 || $"'${head[INHOUSE]}'" ~ /NF/ || NR == 1 ' $FILENAME.alt.NonSyn.exac.gnomad.txt > $FILENAME.alt.NonSyn.exac.gnomad.inhouse.txt




grep -v '^MT' $FILENAME.alt.NonSyn.exac.gnomad.inhouse.txt | grep -v '^GL' | grep -v '^hs' 
> $FILENAME.alt.NonSyn.exac.gnomad.inhouse.noMTnoDecoy.txt


#grep '^MT' $FILENAME.alt.NonSyn.exac.gnomad.inhouse.txt > MT.tmp
#grep '^GL' $FILENAME.alt.NonSyn.exac.gnomad.inhouse.txt > GL.tmp
#grep '^hs' $FILENAME.alt.NonSyn.exac.gnomad.inhouse.txt > hs.tmp
#grep 'CHR:POS:REF:ALT' $FILENAME.alt.NonSyn.exac.gnomad.inhouse.txt > x.tmp
#cat firstline.tmp MT.tmp GL.tmp hs.tmp > $FILENAME.alt.NonSyn.exac.gnomad.inhouse.MT.Decoy.txt
#rm MT.tmp GL.tmp hs.tmp firstline.tmp



# adds a column saying if hom or het
awk -F '\t' -v OFS="\t" '{if ( NR == 1) ajout="Hom"; else if ($0 ~ /1\/1/) ajout="yes"; else ajout="no"; printf ("%s\t%s\n", $0, ajout)}' $FILENAME.alt.NonSyn.exac.gnomad.inhouse.noMTnoDecoy.txt > hom.tmp



awk -F '\t' -v OFS="\t" '{if ( NR == 1) ajout="Het"; else if ($0 ~ /0\/1/ || $0 ~ /1\/0/) ajout="yes"; else ajout="no"; printf ("%s\t%s\n", $0, ajout)}' hom.tmp > $FILENAME.formsc.txt


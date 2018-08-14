#!/bin/bash

#The script produces a file that is used for filtering variants according to the Mutation Significance Cut-off (MSC)  
# and then uses the MSC server to annotate the file

#FILENAME is initialised in all_script
#if you are not using all_scripts then DECOMMENT the next line and write the name of the file
FILENAME="BNvs392.new.merged.vqsr.pass.se"
declare -A header
eval "$(head -n 1 $FILENAME.formsc.txt | awk -F '\t' '{for(i=1;i<=NF;i++) printf "header[%s]=%d ", $i, i }')"

# Sort from the filtered file (*NonSyn.txt) the columns of interest  
awk -F "\t" 'BEGIN{OFS="\t"} {print $"'${header[CHR:POS:REF:ALT]}'"}' $FILENAME.formsc.txt > 1.tmp
awk -F "\t" 'BEGIN{OFS="\t"} {print $"'${header[ID]}'"}' $FILENAME.formsc.txt > 2.tmp
awk -F "\t" 'BEGIN{OFS="\t"} {if ($"'${header[GENE]}'" != "NF") {print $"'${header[GENE]}'"} else {print " "}}' $FILENAME.formsc.txt > 3.tmp
awk -F ":" 'BEGIN{OFS="\t"} {print $1,$2}' 1.tmp > 4.tmp
awk -F ":" 'BEGIN{OFS="\t"} {print $3,$4}' 1.tmp > 5.tmp
paste -d "\t" 4.tmp 2.tmp 5.tmp 3.tmp > 6.tmp
# deletes first line (header)
tail -n +2 6.tmp > $FILENAME.hh.msc

# go and fetch the *.result.txt on 'http://pec630.rockefeller.edu/MSC/UploadServlet'

split -l 3000 $FILENAME.hh.msc $FILENAME.hh.msc.tmp.
for file in $FILENAME.hh.msc.tmp.*
do
	echo "file: $file"
		curl -c cookies.$file.tmp.txt "http://pec630.rockefeller.edu:8080/MSC/UploadServlet"  --referer "http://pec630.rockefeller.edu:8080/MSC/" -H 'Connection: keep-alive' -F "measure1=CADD" -F "measure1=PolyPhen2" -F "measure1=SIFT" -F "confidenceInterval1=ci99" -F "dbSource1=HGMD" -F "displayDataSource1=YES" -F "fname=@./${file}"
		curl -b cookies.$file.tmp.txt "http://pec630.rockefeller.edu:8080/MSC/DownloadServlet"  --referer "http://pec630.rockefeller.edu:8080/MSC/resultVariants.jsp" -H 'Connection: keep-alive' > result.$file.tmp

tail -n +2 result.$file.tmp > result.split.$file.tmp 
done

head -n 1 result.$FILENAME.hh.msc.tmp.aa.tmp > line2.tmp
cat line2.tmp result.split.$FILENAME.hh.msc.tmp.*.tmp > $FILENAME.msc.result.txt.tmp

#sometimes the site doesn't work well and there isn't enough columns
awk -F '\t' -v OFS="\t" '{if ($17 == "") {$17="NF";}}1' $FILENAME.msc.result.txt.tmp > $FILENAME.hh.msc.result.txt

#head is an array containing all the header in order to find the column corresponding to a header (ex: LOF) 
#Nota Bene: we suppose that there are the same annotations for hom and het 
declare -A head
eval "$(head -n 1 $FILENAME.hh.msc.result.txt | awk -F '\t' '{for(i=1;i<=NF;i++) printf "head[%s]=%d ", $i, i }')"

# Sort from the filtered file (*NonSyn.alt.txt) the columns of interest   
awk -F "\t" -v OFS="\t" '{print $6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' $FILENAME.hh.msc.result.txt > 7.tmp

paste -d "\t" $FILENAME.formsc.txt 7.tmp > $FILENAME.hh.alt.NonSyn.exac.gnomad.inhouse.pred.txt.tmp

declare -A title
eval "$(head -n 1 $FILENAME.hh.alt.NonSyn.exac.gnomad.inhouse.pred.txt.tmp | awk -F '\t' '{for(i=1;i<=NF;i++) printf "title[%s]=%d ", $i, i }')"

awk -F '\t' -v OFS='\t' '{if ($"'${title[CLASS]}'" == "TF_binding_site_variant")
{ $"'${title[]}'"=$"'${title[Gene]}'"; $"'${title[GENE_NAME]}'"=$"'${title[Gene]}'";}}1 ' $FILENAME.hh.alt.NonSyn.exac.gnomad.inhouse.pred.txt.tmp > $FILENAME.hh.alt.NonSyn.exac.gnomad.inhouse.pred.txt


mkdir temp
mv *.tmp ./temp
mv *.tmp.* ./temp


#rm *.tmp
#rm *.tmp.*


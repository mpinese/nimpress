#!/bin/bash
#####################################################
## bash script to run preprocess on batch of files ##
#####################################################

## data file needs to be in ../Data/ folder
## second param file read line by line, one for each study, containing relevant info
#DATA="Chang_Hakonarson_2017_Table1.csv"
#PREFIX="Chang_Hakonarson_2017_allparamsv"
## without file extention
#TITLE=${DATA%.*}
#CIT="nDOI: 10.1038/s41467-017-00408-8 PMID: 28924153"
#DEC="Table 1 Association results for the top three genotyped SNPs at 11q22.2"

#echo $TITLE
#cat Data/$DATA
#echo $PREFIX
#echo $TITLE
#echo $CIT

## read in data. What is the easiest? Read lines until marker, for place in comma sperated line and split

while IFS= read -r line || [[ -n "$line" ]]; do
    echo "Text read from file: $line"
    ## split sting
    DATA="$(cut -d',' -f1 <<< $line)"
	PREFIX="$(cut -d',' -f2 <<< $line)"
	## without file extention
	TITLE=${DATA%.*}
	CIT="$(cut -d',' -f3 <<< $line)"
	DEC="$(cut -d',' -f4 <<< $line)"

	echo "Command executed:"
	echo Nimpress_preprocess_pipeline_V5.R ./Data/$DATA --out_prefix $PREFIX --bed /Users/ewilkie/Documents/Polygenic/DataCurationPipeline/mgrb_tier12.bed --title $TITLE --citation $CIT --description $DEC


	Rscript Nimpress_preprocess_pipeline_V5.R ./Data/$DATA --out_prefix "$PREFIX" --bed /Users/ewilkie/Documents/Polygenic/DataCurationPipeline/mgrb_tier12.bed --title "$TITLE" --citation "$CIT" --description "$DEC"

done < $1


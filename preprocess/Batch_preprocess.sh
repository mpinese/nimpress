#!/bin/bash
#####################################################
## bash script to run preprocess on batch of files ##
#####################################################

## data file needs to be in ../Data/ folder
## second param file read line by line, one for each study, containing relevant info
DATA="Chang_Hakonarson_2017_Table1.csv"
PREFIX="Chang_Hakonarson_2017_allparamsv"
## without file extention
TITLE=${DATA%.*}
echo $TITLE

CIT="nDOI: 10.1038/s41467-017-00408-8 PMID: 28924153"


script Nimpress_preprocess_pipeline_V5.R /Data/$DATA --out_prefix $PREFIX --bed /Users/ewilkie/Documents/Polygenic/DataCurationPipeline/mgrb_tier12.bed --title $TITLE --citation $CIT



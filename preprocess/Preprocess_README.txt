
#####################################
## How to run preprocessing script ##
#####################################

#################
## Single file ##
#################

Rscript Nimpress_preprocess_pipeline_V5.R ./Data/$DATA --out_prefix "$PREFIX" --bed /Users/ewilkie/Documents/Polygenic/DataCurationPipeline/mgrb_tier12.bed --title "$TITLE" --citation "$CIT" --description "$DEC"

Where:
$DATA -> .csv file containing rsID, Risk.loci, Freq, OR/BETA, P, Subtype
$PREFIX -> file prefix
$TITLE -> nimpress file title
$CIT -> data publication citation
$DEC -> Data description


##############
## In Batch ##
##############

use the following script
bash Batch_preprocess.sh inparams.txt

where inparams contains one line per file. Each line contains: filename, Prefix, citation and description

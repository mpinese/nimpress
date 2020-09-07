###################################
## Running Nimpress_preprocess.R ##
###################################

Located in the same dir are the following dependencies:
- Nimpress_preprocess_functions.R script containing dependent functions
- Suppl data is located and downloaded to the Suppl folder 
- Output dir Nimpress_preprocess_Output created in this location

Therefore the preprocessing script needs to be invoked in the folder it is downloaded to

General note:
- Only works on genome build version GRCh37
- dbSNP version 151 is used
- tested on R 4.0.0 see details below

###########
## Input ##
###########

Minimum input is a file to be processed in the following format containing header and columns:

rsID,Risk_allele,Freq,<OR> or <BETA>

<rsID> = the dbSNP Reference SNP cluster ID
<Risk_allele> = the allele which has been determined to be associated with the trait/ disease investigated
<Freq> = popultation frequency of the risk allele
<OR> or <BETA> = Odds ratio or beta respectively, numbers will be treated according to header, all value in the header need to be consistent with the type specified (can't mix OR and BETA)

############
## Output ##
############

Finalised data to be used as input to Nimpress is located in the output folder (Nimpress_preprocess_Output) with suffix NIMPRESS_input.txt

Intermediate file that can be used for debugging has suffix Intermediate_results.csv. This file contains the following information:
- dbSNP location, REF and ALT alleles. 
- the GRCh37 allele based on the above location from which strand is inferred
- FLAG.RM indicated whether that rsID has been excluded from the output
- FLAG.AMBIGUOUS indicates whether the REF and ALT alleles are complementary
- BED.COVERAGE whether the rsID is located in the blacklisted region
- INPUT.RISK.TYPE whether risk allele is the REF or the ALT allele
- FLAG.RISK.FLIPPED whether the risk allele is complementary to the REF or ALT and therefore strand flipping has occured
- FLAG.LDPROXY whether a substituted SNV will be used in the output
- LDPROXY.rsID	rsID of substituted SNV
- LDPROXY.CHR	chromosome of substituted SNV
- LDPROXY.START	location of substituted SNV
- LDPROXY.REF.ALLELE	of substituted SNV
- LDPROXY.ALT.ALLELE	of substituted SNV
- INPUT.ALLELE.FREQ	input frequency
- INPUT.BETA input beta (converted from OR if that is supplied instead)

#############################
## Functionality and Usage ##
#############################

The preprocessing script can be run in three different modes. 

The most basic requires only the input file and description and citation parameters 
Rscript Nimpress_preprocess.R --file ./Example/Example_File_to_process.csv --description "Example file" --citation "Authors et al., (2020) Title. Journal"

This will:
- convert OR to Beta via log transformation
- query dbSNP to get the rsID chromosomal locations and reference allele (REF) and alternative allele (ALT)
- rsIDs which don't represent SNVs will be treated as unusable
- BSgenome.Hsapiens.UCSC.hg19 will be used to obtained the strand
- risk allele will be checked against the above to determine whether the allele is ambigious (risk allele is equal to one of REF or ALT and complement of either REF or ALT) and whether strand flipping has occured (risk allele is same as the complement of REF or ALT). Ambiguous alleles are removed from further analysis and won't appear in the output. While strand filling will be corrected.

An extension is enabled by specifying parameters related to a blacklisted bed file containing genomic coordinates of regions that are diffictul to sequence. 

Using the default file (GIAB v2.0 stratification BED file:  https://github.com/genome-in-a-bottle/genome-stratifications)
Rscript Nimpress_preprocess.R --file ./Example/Example_File_to_process.csv --description "Example file" --citation "Authors et al., (2020) Title. Journal" --remove_blacklisted_regions

Using a supplied bed file:
Rscript Nimpress_preprocess.R --file ./Example/Example_File_to_process.csv --description "Example file" --citation "Authors et al., (2020) Title. Journal" --blacklisted_regions_file inhouse_blacklisted_bed.bed

This will:
- perform the same action as basic functionality
- remove any rsIDs form input that are located within the bed regions

Full functionality of the preprocessing script is invoked by providing parameters related to LDproxy:

Rscript Nimpress_preprocess.R --file ./Example/Example_File_to_process.csv --description "Example file" --citation "Authors et al., (2020) Title. Journal" --remove_blacklisted_regions --LDproxy_pop GBR --LDproxy_token abcdef

Rscript Nimpress_preprocess.R --file ./Example/Example_File_to_process.csv --description "Example file" --citation "Authors et al., (2020) Title. Journal" --blacklisted_regions_file inhouse_blacklisted_bed.bed --LDproxy_pop GBR --LDproxy_token abcdef

This will:
- perform the same action as above with the addition of substituting those rsIDs which fall into the blacklisted bed regions with another rsID that has a linkage disequalibrium R-squared value > 0.9
- If no suitable LDproxy substitution rsIDs exist because of R-squared value filtering, all candidates already exist in the input date or are located in the blacklisted bed regions, the rsID will be removed
- SNPs can only be substituted via LDproxy if they are in the 1000 Genome reference panel, a population from that project is selected and a valid API token provided
- monoallelic genes in the population selected will be removed
- 1000 Genomes Project population: https://www.internationalgenome.org/faq/which-populations-are-part-your-study
- generate an LDproxy API token via: https://ldlink.nci.nih.gov/?tab=apiaccess

#########################
## Built and tested on ##
#########################

R version 4.0.0 (2020-04-24)

Libraries:
pacman_0.5.1 
docopt_0.7.1
data.table_1.12.8
GenomicRanges_1.40.0
rentrez_1.2.2
BSgenome.Hsapiens.UCSC.hg19_1.4.3

###############
## Citations ##
###############

Default black-listed genomic regions:
Krusche, P., Trigg, L., Boutros, P.C. et al. 
Best practices for benchmarking germline small-variant calls in human genomes. 
Nat Biotechnol 37, 555–560 (2019). https://doi.org/10.1038/s41587-019-0054-x

LDlinkR:
Myers Timothy A., Chanock Stephen J., Machiela Mitchell J.
LDlinkR: An R Package for Rapidly Calculating Linkage Disequilibrium Statistics in Diverse Populations  
Frontiers in Genetics 11, 157 (2020). https://doi.org/10.3389/fgene.2020.00157   


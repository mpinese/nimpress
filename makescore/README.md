# makescore.R: Convert GWAS summary statistics to nimpress polygenic scores.

## Purpose

`makescore.R` is a script to convert GWAS summary statistics into `nimpress` polygenic score format. Key features of `makescore.R` are:
* Requires only minimal risk allele information as input; suitable for generating a polygenic score from even very basic summary statistics.
* Extensive consistency checking and debug output
* Optional automatic imputation of loci falling into user-supplied blacklisted regions (eg poorly-sequenced regions)

## Requirements

Tested on R version 4.0.3, though it should work on subsequent releases.

This script uses the `renv` library to manage dependencies, for more information see [Installation](#installation).

## Installation

1. If needed, install [R version 4.0.3](https://cran.r-project.org/) or later.
2. Clone the `nimpress` repo.
3. Enter the `makescore` directory.
4. Initialise the environment with the `R` commands:
   ```
   install.packages("renv")
   renv::restore()
   ```

## Usage

`makescore.R` takes as input a table of GWAS summary statistics and produces a polygenic score file that can be read by `nimpress`. `makescore.R` optionally also performs imputation of loci that fall into user-specified blacklisted regions (eg poorly-genotyped regions) by proxy SNPs.

```
makescore.R: Generate an input score file for nimpress from GWAS summary statistics.

Usage:
  makescore.R [options] <infile> <outdir>
  makescore.R (-h | --help)
  makescore.R --version

Arguments:
  <infile>                    Path to the file containing risk loci file in template format (see below).
  <outdir>                    Path to output folder. This will be created if it doesn\'t exist.

General options:
  -h --help                   Show this screen.
  --version                   Show version.

Output options:
  --description="<string>"    Description of the polygenic score, will be copied into the output file. [default: ""]
  --citation="<string>"       Citation of the polygenic score, will be copied into the output file. [default: ""]
  --offset=<offset>           Offset of the polygenic score. This value will be added to calculated scores 
                              (score = sum(dosage*weight) + offset). [default: 0.0]

Blacklist options:
  --blacklist=<file>          BED file of blacklist genome regions. Loci falling in these regions will either be 
                              imputed with proxy SNPs that are not in blacklist regions (if LDproxy_token and 
                              LDproxy_pop are supplied), or dropped (otherwise).

Proxy SNP imputation options:
  --LDproxy_token=<string>    LDproxy API token. Generate at https://ldlink.nci.nih.gov/?tab=apiaccess
  --LDproxy_pop=<string>      Background populations for LDproxy (eg. GBR)
```

### Input file format

The main input is a table that contains summary statistics to transform into a `nimpress` polygenic score. This table is supplied as a tab-separated file with header, with the following fields: `rsID, Risk_allele, Freq, OR or Beta`. Each row corresponds to a locus in the summary statistics. For an example, see `example/wood_height_input_small.tsv`:

rsID      | effectAllele | effectAlleleFrequency | weight
----------|--------------|-----------------------|----------
rs724016  | A            | 0.5551                | -0.0779
rs143384  | A            | 0.5763                | -0.0748
rs8756    | A            | 0.5082                | -0.0585
rs42039   | T            | 0.2655                | 0.0683
rs6845999 | T            | 0.4401                | 0.051


The interpretation of these fields is as follows:
* `rsID`: dbSNP ID for a locus in the summary statistics
* `effectAllele`: the allele which is associated with the trait / disease investigated
* `effectAlleleFrequency`: population frequency of the effect allele
* `weight`: log odds ratio or beta describing the strength of association between the effect allele and the trait investigated.

The `weight` column can contain either log odds ratios or regression coefficients. The final score is calculated as `score = sum(dosage\*weight) + offset`, where `dosage` is the number of copies of the effectAllele in each individual (`dosage \in {0, 1, 2}`). Typically, weight is the reported Beta values (for continuous traits), or log odds (for binary traits / risk), of a GWAS.

### Output files

#### Result output

This script generates two files. An intermediate file for debugging with suffix `.intermediate_results.csv` and a file used as input for NIMPRESS with suffix `.score`. The input file name will be used as prefix. The example file will therefore generate the following output files:

#### More details on intermetidate results

The intermediate file has the following columns that will aid in understanding the changes that have been applied to the input data to generate the output data

- dbSNP location, REF and ALT alleles. 
- the GRCh37 allele based on the above location from which strand is inferred
- FLAG.RM indicated whether that rsID has been excluded from the output
- FLAG.AMBIGUOUS indicates whether the REF and ALT alleles are complementary
- BED.COVERAGE whether the rsID is located in the blacklisted region
- INPUT.RISK.TYPE whether risk allele is the REF or the ALT allele
- FLAG.RISK.FLIPPED whether the risk allele is complementary to the REF or ALT and therefore strand flipping has occured
- FLAG.LDPROXY whether a substituted SNV will be used in the output
- LDPROXY.rsID  rsID of substituted SNV
- LDPROXY.CHR chromosome of substituted SNV
- LDPROXY.START location of substituted SNV
- LDPROXY.REF.ALLELE  of substituted SNV
- LDPROXY.ALT.ALLELE  of substituted SNV
- INPUT.ALLELE.FREQ input frequency
- INPUT.BETA input beta (converted from OR if that is supplied instead)

## Examples

### 1. Basic
In the most basic application, `makescore.R` requires only the GWAS summary statistics:
`Rscript makescore.R example/wood_height_input_small.tsv example/wood_height_small_example1`

This will:
- convert OR to Beta via log transformation
- query dbSNP to get the rsID chromosomal locations and reference allele (REF) and alternative allele (ALT)
- rsIDs which don't represent SNVs will be treated as unusable
- BSgenome.Hsapiens.UCSC.hg19 will be used to obtain the strand
- effect allele will be checked against the above to determine whether the allele is ambigious (effect allele is equal to one of REF or ALT and complement of either REF or ALT) and whether strand flipping has occured (effect allele is same as the complement of REF or ALT). Ambiguous alleles are removed from further analysis and won't appear in the output. Unambiguously strand flipped loci will be corrected.

### 2. With GBS blacklist

An extension is enabled by specifying parameters related to a blacklisted bed file containing genomic coordinates of regions that are difficult to sequence. 

`Rscript makescore.R --blacklist=example/wood_height_blacklist.bed example/wood_height_input_small.tsv example/wood_height_small_example2`

This will:
- perform the same action as basic functionality
- remove any rsIDs from the input that are located within the blacklist bed regions

### 3. With GBS blacklist and imputation

Full functionality of the preprocessing script is invoked by providing parameters related to LDproxy:

`Rscript makescore.R --blacklist=example/wood_height_blacklist.bed --LDproxy_token=<token> --LDproxy_pop=GBR example/wood_height_input_small.tsv example/wood_height_small_example3`

This will:
- perform the same action as above with the addition of substituting those rsIDs which fall into the blacklisted bed regions with another rsID that has a linkage disequalibrium R-squared value > 0.9
- If no suitable LDproxy substitution rsIDs exist because of R-squared value filtering, all candidates already exist in the input date or are located in the blacklisted bed regions, the rsID will be removed
- SNPs can only be substituted via LDproxy if they are in the 1000 Genome reference panel, a population from that project is selected, and a valid API token provided
- monoallelic genes in the population selected will be removed
- 1000 Genomes Project population: https://www.internationalgenome.org/faq/which-populations-are-part-your-study
- generate an LDproxy API token via: https://ldlink.nci.nih.gov/?tab=apiaccess

## Notes and Limitations

* Currently this script is fixed on dbSNP v151, and genome build version GRCh37
* Internet access is required, and the script is quite slow (~ 3 seconds / locus) due to web access to dbSNP. A future version will use a local cache of dbSNP for speed.
* Running `Rscript makescore.R` gives a truncated output; this is a limitation of the `docopt` implementation in R. To get a full command line help message, run `Rscript makescore.R --help`.

## References

LDlinkR:
Myers Timothy A., Chanock Stephen J., Machiela Mitchell J.
LDlinkR: An R Package for Rapidly Calculating Linkage Disequilibrium Statistics in Diverse Populations  
Frontiers in Genetics 11, 157 (2020). https://doi.org/10.3389/fgene.2020.00157   


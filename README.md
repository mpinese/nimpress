# nimpress: Polygenic score calculation direct from VCF/BCF

```
  Usage:
    nimpress [options] <scoredef> <genotypes.vcf>
    nimpress (-h | --help)
    nimpress --version

  Options:
    -h --help         Show this screen.
    --version         Show version.
    --cov=<path>      Path to a BED file supplying genome regions that have been
                      genotyped in the genotypes.vcf file.
    --imp-locus=<m>   Imputation to apply for whole loci which are either not
                      in the sequenced BED regions, or fail (QUAL flag or too
                      many samples with missing genotype). Valid values are ps, 
                      homref, fail [default: ps].
    --imp-sample=<m>  Imputation to apply for an individual sample with missing 
                      genotype. Valid values are ps, homref, fail, int_fail, 
                      int_ps [default: int_ps].
    --maxmis=<f>      Maximum fraction of samples with missing genotypes allowed
                      at a locus.  Loci containing more than this fraction of 
                      samples missing will be considered bad, and have all 
                      genotypes (even non-missing ones) imputed [default: 0.05].
    --mincs=<n>       Minimum number of genotypes.vcf samples without missing 
                      genotype at a locus for this locus to be eligible for 
                      internal imputation [default: 100].
    --afmisp=<f>      p-value threshold for warning about allele frequency 
                      mismatch between the polygenic score and the supplied 
                      cohort [default: 0.001].

  Imputation methods:
  ps        Impute with dosage based on the polygenic score effect allele 
            frequency.
  homref    Impute to homozygous reference genotype.
  fail      Do not impute, but fail. Failed samples will have a score of "nan"
  int_ps    Impute with dosage calculated from non-missing samples in the 
            cohort. At least --mincs non-missing samples must be available for 
            this method to be used, else it will fall back to ps.
  int_fail  Impute with dosage calculated from non-missing samples in the 
            cohort. At least --mincs non-missing samples must be available for 
            this method to be used, else it will fall back to fail.
```

# Polygenic score format
At the moment a placeholder format is used for polygenic scores. This is likely to change rapidly.

The placeholder format is a text file:
```
<name>
<description>
<citation>
<genome version>
<offset>
<chrom>\t<pos>\t<ref>\t<alt>\t<beta>\taaf
<chrom>\t<pos>\t<ref>\t<alt>\t<beta>\taaf
...
```
The first four header records (name, description, citation, and genome version) are free text, 
the last header record (offset) is a string representation of a floating point number.

Lines following the header define the polygenic score alleles and coefficients, one line per
allele, with fields separated by tabs.  <beta> is the polygenic score coefficient associated with
a single alternate allele.  <aaf> is the alt allele frequency in the polygenic score derivation
cohort, 0 < aaf <= 1.

The polygenic score for sample `i` is calculated as:
```
  score_i = (offset + sum_{j=1..m}(beta_j*dosage_ij)) / m
```
where `offset` is the offset as given in the header, `beta_j` is the beta for row `j` of the 
polygenic score definition, `dosage_ij` is the dosage of the row `j` alt allele in sample `i`,
and `m` is the number of alleles (rows) in the polygenic score.  `dosage_ij` may be imputed.

# Limitations
* Currently the coverage BED file is ignored.
* Allele frequency mismatch tests are not yet implemented.
* Currently has diploid-specific logic.

# Future
* Option to calculate soft PS based on genotype likelihoods (VCF PL field)

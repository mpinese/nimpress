# nimpress: Polygenic score calculation direct from VCF/BCF

[![Build Status](https://travis-ci.org/mpinese/nimpress.svg?branch=master)](https://travis-ci.org/mpinese/nimpress)

**nimpress** is a lightweight tool to calculate polygenic scores from genomic data. It has a simple score input format (and comes with a utility to generate scores from published GWAS summary statistics), works directly from VCF and BCF with no intermediate files, and has no special requirements beyond the *htslib* library.

## Prerequisites

**nimpress** is supplied as a pre-built static-linked binary which runs on modern Linux systems with no requirements. **nimpress** is also available as a Docker container that can run on Windows, Mac, and Unix/Linux systems with the Docker runtime installed.

The static binary or Docker container will suit most users. However, if you wish to build **nimpress** from source, ensure you have a working installation of the *nim* language (https://nim-lang.org/), version >= 1.0.0, with packages:
* `docopt` >= 0.6.8
* `hts` >= 0.2.21
* `lapper` >= 0.1.5

A working installation of *htslib* (http://www.htslib.org/) is also required.

The easiest way to set up a suitable *nim* environment is using *nimble*:
1. Ensure you have a working installation of *htslib* (http://www.htslib.org/).
2. Install *nim* following the directions at https://nim-lang.org/. This will also install the *nimble* package manager
3. Clone this repository `git clone https://github.com/mpinese/nimpress.git`
4. Inside the **nimpress** directory run `nimble install -d` to install all dependencies.

## Installation

### Static pre-built binary (Linux only)
1. Download the latest pre-built binary from the [releases page](https://github.com/mpinese/nimpress/releases "Nimpress Releases Page"), and save to a location of your choice.
2. Depending on your system, you may need to mark the downloaded file as executable with `chmod +x nimpress-x86_64`
3. Verify a successful installation with `./nimpress-x86-64`; you should see a brief usage message.

### Docker container (Linux, Mac, Windows)
1. Install the Docker runtime and open a Docker console. See https://www.docker.com/ for instructions on how to do this for each system.
2. Within the Docker console, fetch the **nimpress** docker image with `docker pull mpinese/nimpress`
3. Verify a successful installation with `docker run -t mpinese/nimpress`; you should see a brief usage message.

### From source (Linux, Mac, Windows)
1. Install *nim* following the instructions in [Prerequisites](#prerequisites).
2. Clone this repository: `git clone https://github.com/mpinese/nimpress.git`
3. From the repository directory (`cd nimpress`) run `nimble install`.
4. Verify installation with `nimpress`; you should see a brief usage message.

## Usage

Specifics of usage depend on how **nimpress** was installed, either as a native (source install or static binary) or Docker install:

### Source install or static binary
**nimpress** is run as a standard binary from the command line, with `nimpress` or `./nimpress-x86_64`. Command line options are supplied as detailed in [Shared options](#shared-options) below.

### Docker
If running **nimpress** inside a Docker container, you must supply input files (genotype vcf/bcf with its index, and polygenic score file) to the container so that **nimpress** can access them. This can be accomplished in the `docker run` command as:
```
docker run \
  -v <local_path_to_score_file>:/root/ps.scores \
  -v <local_path_to_genotype_file>:/root/genotypes.<bcf|vcf.gz> \
  -v <local_path_to_genotype_index>:/root/genotypes.<bcf|vcf.gz>.<csi|tbi> \
  -t mpinese/nimpress <nimpress options here> /root/ps.scores /root/genotypes.<bcf|vcf.gz>
```
A coverage BED, if used, would also be supplied with a `-v` flag. Other nimpress options would be supplied following the [Shared options](#shared-options) section below.

For example, if you have genotypes in a bcf file `/home/user/mycohort.bcf` with index `/home/user/mycohort.bcf.csi`, a genotyping coverage bed in `/home/user/genotyped_regions.bed`, and wish to assess the Wood *et al* height PS in `/home/user/nimpress/scores/wood-25282103-height.scores`, you would use the command:
```
docker run \
  -v /home/user/nimpress/scores/wood-25282103-height.scores:/root/ps.scores \
  -v /home/user/mycohort.bcf:/root/genotypes.bcf \
  -v /home/user/mycohort.bcf.csi:/root/genotypes.bcf.csi \
  -v /home/user/genotyped_regions.bed:/root/genotyped.bed \
  -t mpinese/nimpress --cov=/root/genotyped.bed /root/ps.scores /root/genotypes.bcf
```
Note that if running on Windows, full paths (eg `C:\Users\User\Desktop\nimpress\scores\wood-25282103-height.scores`) must be supplied, or else an error like `string=is not a valid Windows path` will be encountered.

### Shared options
For either installation method, **nimpress** is configured with the following command line flags:
```
  Usage:
    nimpress [options] <scoredef> <genotypes.bcf|genotypes.vcf.gz>
    nimpress (-h | --help)
    nimpress --version

  Options:
    -h --help          Show this screen.
    --version          Show version.
    --cov=<path>       Path to a BED file supplying genome regions 
                       that have been genotyped in the genotypes file.
    --imp-locus=<m>    Imputation to apply for whole loci which are 
                       either not in the sequenced BED regions, or 
                       fail (too many samples with missing genotype, 
                       as set by --maxmis, or failing VCF QUAL field 
                       if --ignorefilt is not set). Valid values are 
                       ps, homref, fail, ignore [default: ps].
    --imp-missing=<m>  Imputation to apply for loci which are in the 
                       sequenced BED regions (and thus should have 
                       been genotyped), but are completely missing 
                       from the VCF. Valid values are homref, ignore 
                       [default: homref].
    --imp-sample=<m>   Imputation to apply for an individual sample 
                       with missing genotype. Valid values are ps, 
                       homref, fail, int_fail, int_ps 
                       [default: int_ps].
    --maxmis=<f>       Maximum fraction of samples with missing 
                       genotypes allowed at a locus. Loci containing 
                       more than this fraction of samples missing will
                       be considered bad, and have all genotypes (even
                       non-missing ones) imputed [default: 0.05].
    --mincs=<n>        Minimum number of genotypes.vcf samples without
                       missing genotype at a locus for this locus to 
                       be eligible for internal imputation [default: 
                       100].
    --afmisp=<f>       p-value threshold for warning about allele 
                       frequency mismatch between the polygenic score 
                       and the supplied cohort [default: 0.001].
    --ignorefilt       Ignore the VCF FILTER field. If set, all 
                       variants in the VCF will be used regardless of 
                       FILTER field contents. If not set, variants 
                       with a FILTER field other than "." or "PASS" 
                       will always be imputed.

  Imputation methods:
  ps        Impute with dosage based on the polygenic score effect 
            allele frequency.
  homref    Impute to homozygous reference genotype.
  fail      Do not impute, but fail. Failed samples will have a score 
            of "nan".
  ignore    Completely ignore missing loci, as if they were never in 
            the score definition.
  int_ps    Impute with dosage calculated from non-missing samples in 
            the cohort. At least --mincs non-missing samples must be 
            available for this method to be used, else it will fall 
            back to ps.
  int_fail  Impute with dosage calculated from non-missing samples in 
            the cohort. At least --mincs non-missing samples must be 
            available for this method to be used, else it will fall 
            back to fail.
```

# Polygenic score format
Polygenic scores are described by text files with the following format:
```
<name>
<description>
<citation>
<genome version>
<offset>
<chrom>\t<pos>\t<refallele>\t<effallele>\t<beta>\t<effallele_af>
<chrom>\t<pos>\t<refallele>\t<effallele>\t<beta>\t<effallele_af>
...
```
The first four header records (name, description, citation, and genome version) are free text, the last header record (offset) is a string representation of a floating point number.

Lines following the header define the polygenic score alleles and coefficients, one line per allele, with fields separated by tabs. <beta> is the polygenic score coefficient associated with a single alternate allele. <effallele_af> is the effect allele frequency in the polygenic score derivation cohort, 0 < effallele_af <= 1.

The polygenic score for sample `i` is calculated as:
```
  score_i = sum_{j=1..m}(beta_j*dosage_ij)/m + offset
```
where `offset` is the offset as given in the header, `beta_j` is the beta for row `j` of the polygenic score definition, `dosage_ij` is the dosage of the row `j` effect allele in sample `i`, and `m` is the number of alleles (rows) in the polygenic score.  `dosage_ij` may be imputed. It's not uncommon for the reference allele to also be the effect allele; in this case set `<effallele>` to equal `<refallele>`.

# Limitations
* Currently diploid-specific.
* Does not fully handle multi-allelic risk loci (specifically, loci at which more than one allele has a nonzero beta are not supported)
* Performs only simple allele matching. As the representation of some variants in VCF is not unique, this may lead to polygenic score variants being imputed even if they are present in the VCF. A partial workaround is to pass all variants (both in the genotype file and in the polygenic score file) through variant normalization (using `vt norm` or similar, https://genome.sph.umich.edu/wiki/Vt) before calculating polygenic scores, this is effective in most cases.

# Future
* Option to calculate soft PS based on genotype likelihoods (VCF PL field)
* Full multi-allelic locus handling (will require imputation to a distribution over dosages of all alleles, instead of just the effect allele)
* Smarter variant matching.

# Contributing
Contributions are encouraged. Please fork this repository, make changes in your forked version, and submit a pull request.

# Acknowledgements
* [Emilie Wilkie](https://github.com/ewilkie/) for preprocessing code & much debugging.
* [Brent Pedersen](https://github.com/brentp/) for the excellent hts-nim and lapper libraries.

# Contact
Mark Pinese, Computational Biology Group, Children's Cancer Institute (<MPinese@ccia.org.au>)

# Licence
MIT

import logging
import macros
import math
import strUtils
import tables

import docopt
import hts
import lapper



################################################################################
## Utility functions
################################################################################


macro tpub*(x: untyped): untyped =
  ## Marks a proc with an export asterisk when ``-d:testing`` is defined.
  ## From https://forum.nim-lang.org/t/3269
  expectKind(x, RoutineNodes)
  when defined(testing):
    let n = name(x)
    x.name = newTree(nnkPostfix, ident"*", n)
  result = x


proc isNaN(x: float): bool =
  result = x.classify == fcNaN


proc tallyAlleles(rawDosages: seq[float]): (float, float, float) =
  ## Tally the alleles in rawDosages, as returned by getRawDosages.  Returns
  ## a tuple of three floats, with entries:
  ##   number of samples with genotype
  ##   number of samples missing genotype
  ##   total count of effect allele in genotyped samples.
  var ngenotyped = 0.0
  var nmissing = 0.0
  var neffectallele = 0.0
  for dosage in rawDosages:
    if dosage.isNaN:
      nmissing += 1.0
    else:
      ngenotyped += 1.0
      neffectallele += dosage
  return (ngenotyped, nmissing, neffectallele)


proc dbinom (x: int, n: int, p: float): float {.tpub.} =
  ## Pr(x successes in n trials | Pr(success) = p)
  binom(n, x).toFloat * pow(p, x.toFloat) * pow(1.0-p, (n-x).toFloat)


proc pbinom(x: int, n: int, p: float): float {.tpub.} =
  ## Returns the lower tail binomial probability of x successes or fewer in n
  ## trials, given a success probability of p: Pr(k <= x; n, p)
  result = 0.0
  for xi in 0..x:
    result += dbinom(xi, n, p)


proc binomTest(x: int, n: int, p: float): float {.tpub.} =
  ## Two-sided binomial test of observing x successes or more extreme in n
  ## trials, given success probability of p.  Returns the p value.

  # Edge cases
  #if p == 0.0:
  #  return if x == 0: 1.0 else: 0.0
  #elif p == 1.0:
  #  return if x == n: 1.0 else: 0.0

  #let
  #  probx = dbinom(x, n, p)
  #  expectedVal = n.toFloat*p

  # Edge case again
  #if abs(x.toFloat/expectedVal - 1.0) < 1.0e-6:
  #  return 1.0

  # Find the integration limits by enumeration, and use to evaluate the
  # p value. There's probably a more efficient implementation (Newton's?) for
  # large n; consider this if profiling indicates this is a bottleneck.
  #if x.toFloat < expectedVal:
  #  var y = 0
  #  for xi in expectedVal.ceil.toInt..n:
  #    if dbinom(xi, n, p) <= probx * (1.0+1.0e-7):
  #      y += 1
  #  return pbinom(x, n, p) + (1.0 - pbinom(n - y, n, p))
  #else: # x > expectedVal
  #  var y = 0
  #  for xi in 0..(floor(expectedVal).toInt):
  #    if dbinom(xi, n, p) <= probx * (1.0+1.0e-7):
  #      y += 1
  #  return pbinom(y - 1, n, p) + (1.0 - pbinom(x - 1, n, p))
  return 1.0


################################################################################
## Polygenic score file object
################################################################################

type ScoreFile* = object
  ## Really rough polygenic score file definition, just to get something
  ## working.  Current format is 5 header lines followed by 6-column TSV.
  ## Header lines are:
  ##   name (string)
  ##   description (string)
  ##   citation (string)
  ##   genome version (string)
  ##   offset (string representation of float)
  ## The subsequent TSV section is headerless, with one row per effect allele,
  ## columns:
  ##   chrom, pos, ref, effectallele, beta, eaf
  ## where beta is the PS coefficient and eaf the effectallele MAF in the source
  ## population.  effectallele may equal ref, in which case beta is the
  ## coefficient for reference allele dosage.
  ##
  ## In future this should be a 'real' format (tabix-compatible? Will need to be
  ## space efficient if genome-wide scores are on the table).
  fileobj: File
  name: string
  desc: string
  cite: string
  genomever: string
  offset: float


type ScoreEntry = tuple
  ## Tuple container for ScoreFile records
  contig: string
  pos: int
  refseq: string
  easeq: string
  beta: float
  eaf: float

proc stop*(s:ScoreEntry): int {.inline.} =
  return s.pos + s.refseq.len - 1

proc open*(scoreFile: var ScoreFile, path: string): bool =
  ## Open a ScoreFile
  scoreFile.fileobj = open(path)
  if scoreFile.fileobj.isNil:
    return false
  scoreFile.name = scoreFile.fileobj.readLine.strip(leading = false)
  scoreFile.desc = scoreFile.fileobj.readLine.strip(leading = false)
  scoreFile.cite = scoreFile.fileobj.readLine.strip(leading = false)
  scoreFile.genomever = scoreFile.fileobj.readLine.strip(leading = false)
  scoreFile.offset = scoreFile.fileobj.readLine.strip(
      leading = false).parseFloat
  return true


iterator items(scoreFile: ScoreFile): ScoreEntry =
  ## Iterate over entries in scoreFile
  var line: string
  while scoreFile.fileobj.readLine(line):
    let lineparts = line.strip(leading = false).split('\t')
    doAssert lineparts.len == 6
    yield (lineparts[0], lineparts[1].parseInt, lineparts[2], lineparts[3],
           lineparts[4].parseFloat, lineparts[5].parseFloat)



################################################################################
## Handling of well-genotyped regions
################################################################################

type ContigInterval = object
  start: int
  stop: int

proc start(ival: ContigInterval): int = return ival.start
proc stop(ival: ContigInterval): int = return ival.stop
proc `$`(ival: ContigInterval): string = return "[$#,$#)" %
  [$ival.start, $ival.stop]


type GenomeIntervals* = object
  init: bool
  contigIntervalsIndex: Table[string, Lapper[ContigInterval]]
  contigIntervals: Table[string, seq[ContigInterval]]


proc loadBedIntervals*(ivals: var GenomeIntervals, path: string): bool =
  ## Load a BED file in path into the ivals object. On success, all ivals data
  ## will be replaced with the BED file intervals.  Returns true on success,
  ## else false.  If the load fails (return value false), ivals will not be
  ## changed.
  let infile = open(path)
  if infile.isNil:
    return false

  # Clear any existing intervals, and load the BED ones
  ivals.init = false
  ivals.contigIntervalsIndex.clear()
  ivals.contigIntervals.clear()
  for line in infile.lines:
    let lineparts = line.strip(leading = false).split('\t')
    doAssert lineparts.len >= 3
    let
      chrom = lineparts[0]
      start0 = lineparts[1]
      end1 = lineparts[2]
    if not ivals.contigIntervals.hasKey(chrom):
      ivals.contigIntervals[chrom] = @[]
    ivals.contigIntervals[chrom].add(ContigInterval(start: start0.parseInt,
        stop: end1.parseInt))

  # Add lapper indices
  for chrom in ivals.contigIntervals.keys:
    ivals.contigIntervalsIndex[chrom] = lapify(ivals.contigIntervals[chrom])

  ivals.init = true
  return true

proc contains(ci:ContigInterval, scoreEntry: ScoreEntry): bool {.inline.} =
  return ci.start < scoreEntry.pos and ci.stop >= scoreEntry.stop

proc isVariantCovered(scoreEntry: ScoreEntry,
    coveredIvals: GenomeIntervals): bool =
  ## Tests for whether the variant in scoreEntry falls within an interval in
  ## coveredIvals.  coveredIvals contains half-open 0-based intervals, as per
  ## the BED format.  scoreEntry.pos is in 1-based format, as per VCF. If
  ## coveredIvals is not initalised (coveredIvals.init == false), always returns
  ## false.
  ##
  ## Return value:
  ##   true if scoreEntry falls entirely within at least one coveredIvals
  ##     interval,
  ##   else false.
  if not coveredIvals.contigIntervalsIndex.hasKey(scoreEntry.contig):
    log(lvlWarn, "Contig " & scoreEntry.contig & " not present within the " &
        "coverage BED file.")
    return false # The whole contig is missing from coveredIvals.

  # Get all overlapping intervals. Note the -1 to convert from 1-based inclusive
  # VCF pos to BED half-open intervals
  if scoreEntry.contig notin coveredIvals.contigIntervalsIndex:
    return false

  var ivs = new_seq[ContigInterval]()
  var contigLapper = coveredIvals.contigIntervalsIndex[scoreEntry.contig]
  let anyfound = contigLapper.find(scoreEntry.pos-1, scoreEntry.stop + 1, ivs)
  if not anyfound:
    return false

  # At least one overlapping interval was found. Loop over each to confirm the
  # variant falls entirely within the interval.
  for ci in ivs:
    if ci.contains(scoreEntry): return true
  return false



################################################################################
## VCF access: variant search and dosage querying
################################################################################

proc findVariant(contig: string, pos: int, refseq: string, easeq: string,
                 vcf: VCF): Variant =
  # Find contig:pos:refseq:easeq in vcf.  Returns the
  # whole VCF Variant if found, else nil.
  result = nil
  for variant in vcf.query(contig & ":" & $pos & "-" & $(pos + refseq.len - 1)):
    if variant.REF == refseq:
      if easeq == refseq:
        return variant
      for valt in variant.ALT:
        if valt == easeq:
          return variant


proc getRawDosages(rawDosages: var seq[float], variant: Variant,
    easeq: string) =
  ## Get the dosages of easeq in the VCF Variant variant.  Returns
  ## a sequence with values in {NaN, 0., 1., 2.}, being the dosage, or NaN
  ## if no genotype is available.
  ## TODO: Had trouble figuring out the best way to access the hts-nim API for
  ## this.  Probably a much faster way to do it.  Not recreating the gts int32
  ## seq each time seems like a good start.
  let eaidx =
    if easeq == variant.REF:
      0
    else:
      find(variant.ALT, easeq) + 1 # index of the desired alt allele
  doAssert eaidx >= 0
  var gts = newSeqUninitialized[int32](variant.n_samples)

  var i = 0
  for gt in genotypes(variant.format, gts):
    rawDosages[i] = 0.0
    for allele in gt:
      if value(allele) == eaidx:
        rawDosages[i] += 1
      elif value(allele) == -1:
        rawDosages[i] = NaN
    i += 1



################################################################################
## Imputation
################################################################################

## Locus and sample imputation methods:
## ps        Impute with dosage based on the polygenic score effect allele
##           frequency.
## homref    Impute to homozygous reference genotype.
## fail      Do not impute, but fail. Failed samples will have a score of "nan"
## ignore    Completely ignore missing loci, as if they were never in the score 
##           definition.
## int_ps    Impute with dosage calculated from non-missing samples in the
##           cohort. At least --mincs non-missing samples must be available for
##           this method to be used, else it will fall back to ps.
## int_fail  Impute with dosage calculated from non-missing samples in the
##           cohort. At least --mincs non-missing samples must be available for
##           this method to be used, else it will fall back to fail.
type ImputeMethodLocus*{.pure.} = enum ps, homref, fail, ignore
type ImputeMethodMissing*{.pure.} = enum homref, ignore
type ImputeMethodSample*{.pure.} = enum ps, homref, fail, int_ps, int_fail


proc imputeLocusDosages(dosages: var seq[float], scoreEntry: ScoreEntry,
                        imputeMethodLocus: ImputeMethodLocus): bool =
  ## Impute all dosages at a locus.  Even non-missing genotypes are imputed.
  ##
  ## dosages: destination seq to which imputed dosages will be written.
  ## scoreEntry: the polygenic score entry corresponding to this locus.
  ## imputeMethodLocus: imputation method.
  ##
  ## Returns a boolean indicating if the locus should be used for PRS 
  ## calculation (return value true), or if it should be discarded (return value
  ## false).
  if imputeMethodLocus == ImputeMethodLocus.ignore:
    return false

  let imputed_dosage = case imputeMethodLocus:
    of ImputeMethodLocus.ps:
      scoreEntry.eaf*2.0
    of ImputeMethodLocus.homref:
      if scoreEntry.refseq == scoreEntry.easeq:
        2.0
      else:
        0.0
    of ImputeMethodLocus.fail:
      NaN
    of ImputeMethodLocus.ignore:
      NaN

  for i in 0..dosages.high:
    dosages[i] = imputed_dosage

  return true


proc imputeSampleDosages(dosages: var seq[float], scoreEntry: ScoreEntry,
                         nEffectAllele: float, nGenotyped: float,
                         minGtForInternalImput: int,
                         imputeMethodSample: ImputeMethodSample) =
  ## Impute missing genotype dosages at a locus.  Genotypes which are not
  ## missing will not be imputed.
  ##
  ## dosages: destination seq to which imputed dosages will be written.
  ## scoreEntry: the polygenic score entry corresponding to this locus.
  ## imputeMethodSample: imputation method.
  let imputed_dosage = case imputeMethodSample:
    of ImputeMethodSample.ps:
      scoreEntry.eaf*2.0
    of ImputeMethodSample.homref:
      if scoreEntry.refseq == scoreEntry.easeq:
        2.0
      else:
        0.0
    of ImputeMethodSample.fail:
      NaN
    of ImputeMethodSample.int_ps, ImputeMethodSample.int_fail:
      if nGenotyped >= minGtForInternalImput.toFloat:
        nEffectAllele / nGenotyped
      else:
        if imputeMethodSample == ImputeMethodSample.int_ps:
          scoreEntry.eaf*2.0
        else:
          NaN

  for i in 0..dosages.high:
    if dosages[i].isNaN:
      dosages[i] = imputed_dosage


proc getImputedDosages(dosages: var seq[float], scoreEntry: ScoreEntry,
                       genotypeVcf: VCF, restrictToCoveredRgns: bool,
                       coveredIvals: GenomeIntervals,
                       imputeMethodLocus: ImputeMethodLocus,
                       imputeMethodMissing: ImputeMethodMissing,
                       imputeMethodSample: ImputeMethodSample,
                       maxMissingRate: float, afMismatchPthresh: float,
                       minGtForInternalImput: int, 
                       ignoreFilterField: bool): bool =
  ## Fetch dosages of allele described in scoreEntry from samples genotyped in
  ## genotypeVcf.  Impute dosages if necessary.
  ##
  ##   dosages: destination seq to which imputed dosages will be written.
  ##   scoreEntry: the polygenic score entry corresponding to this locus.
  ##   restrictToCoveredRgns: true if analysis should be restricted to
  ##                          well-called regions in coveredIvals; false to
  ##                          consider all loci to be well-covered.
  ##   coveredIvals: object containing genome regions which have been well-
  ##                 called (covered) by the genotyping method.
  ##   imputeMethodLocus: Imputation method to use when a whole locus fails or
  ##                      is missing / not covered.
  ##   imputeMethodMissing: Imputation to use for variants which are at covered 
  ##                        loci, but which are not present in the VCF.
  ##   imputeMethodSample: Imputation method to use for individual samples with
  ##                       missing genotype, in a locus that passes QC filters.
  ##   maxMissingRate: loci with more than this rate of missing samples fail
  ##                   QC and are imputed.
  ##   afMismatchPthresh: p-value threshold to warn about allele frequency
  ##                      mismatch between the cohort in genotypeVcf and the
  ##                      polygenic score in scoreFile.
  ##   minGtForInternalImput: Minimum number of genotyped samples at a locus for
  ##                          internal imputation to be applied.
  ##   ignoreFilterField: If true, do not consider the VCF FILTER field when
  ##                      identifying loci for imputation. If false, loci with
  ##                      FILTER other than "." or "PASS" will always be imputed.
  ##
  ## Returns a boolean indicating if the locus should be used for PRS 
  ## calculation (return value true), or if it should be discarded (return value
  ## false).
  let nsamples = genotypeVcf.n_samples
  dosages.setLen(nsamples)

  if restrictToCoveredRgns and not isVariantCovered(scoreEntry, coveredIvals):
    log(lvlWarn, "Locus " & scoreEntry.contig & ":" & $scoreEntry.pos & "-" &
        $scoreEntry.stop &
        " is not covered by the sequence coverage BED.  Imputing all dosages " &
        "at this locus.")
    return imputeLocusDosages(dosages, scoreEntry, imputeMethodLocus)
  
  let variant = findVariant(scoreEntry.contig, scoreEntry.pos,
                            scoreEntry.refseq, scoreEntry.easeq, genotypeVcf)
  
  if variant.isNil:
    if binomTest(0, nsamples*2, scoreEntry.eaf) < afMismatchPthresh:
      log(lvlWarn, "Variant " & scoreEntry.contig & ":" & $scoreEntry.pos &
          ":" & $scoreEntry.refseq & ":" & $scoreEntry.easeq &
          " cohort EAF is 0 in " & $nsamples & " samples.  This is highly" &
          " unlikely given polygenic score EAF of " & $scoreEntry.eaf)
    # This variant is in the covered regions (or no coverage BED was supplied,
    # in which case we assume it's covered), but is missing from the VCF. 
    # Impute as specified by imputeMethodMissing.
    if imputeMethodMissing == ImputeMethodMissing.homref:
      let impute_dosage = if scoreEntry.refseq == scoreEntry.easeq: 2.0 else: 0.0
      for i in 0..dosages.high:
        dosages[i] = impute_dosage
      return true
    elif imputeMethodMissing == ImputeMethodMissing.ignore:
      return false
  
  if ignoreFilterField == false and $variant.FILTER != "." and $variant.FILTER != "PASS":
    log(lvlWarn, "Variant " & scoreEntry.contig & ":" & $scoreEntry.pos & ":" &
        $scoreEntry.refseq & ":" & $scoreEntry.easeq &
        " has a FILTER flag set (value \"" & $variant.FILTER & "\").  " &
        "Imputing all dosages at this locus.")
    return imputeLocusDosages(dosages, scoreEntry, imputeMethodLocus)

  # Fetch the raw dosages (values in {NaN, 0., 1., 2.}) from the VCF.
  getRawDosages(dosages, variant, scoreEntry.easeq)

  let (ngenotyped, nmissing, neffectallele) = tallyAlleles(dosages)

  let missingrate = nmissing / nsamples.toFloat
  if missingrate > maxMissingRate:
    log(lvlWarn, "Locus " & scoreEntry.contig & ":" & $scoreEntry.pos & "-" &
        $scoreEntry.stop & " has " &
        $(missingrate*100) & "% of samples missing a genotype. This exceeds " &
        "the missingness threshold; imputing all dosages at this locus.")
    return imputeLocusDosages(dosages, scoreEntry, imputeMethodLocus)

  if binomTest(neffectallele.toInt, (nsamples-nmissing.toInt)*2,
      scoreEntry.eaf) < afMismatchPthresh:
    log(lvlWarn, "Variant " & scoreEntry.contig & ":" & $scoreEntry.pos &
        ":" & $scoreEntry.refseq & ":" & $scoreEntry.easeq & " cohort EAF is " &
        $(neffectallele/((nsamples-nmissing.toInt)*2).toFloat) & " in " &
        $nsamples & " samples.  This is highly unlikely given polygenic " &
        "score EAF of " & $scoreEntry.eaf)

  # Impute single missing sample dosages
  imputeSampleDosages(dosages, scoreEntry, neffectallele, ngenotyped,
                      minGtForInternalImput, imputeMethodSample)

  return true


################################################################################
## Polygenic score calculation
################################################################################

proc computePolygenicScores*(scores: var seq[float], scoreFile: ScoreFile,
                             genotypeVcf: VCF, restrictToCoveredRgns: bool,
                             coveredIvals: GenomeIntervals,
                             imputeMethodLocus: ImputeMethodLocus,
                             imputeMethodMissing: ImputeMethodMissing,
                             imputeMethodSample: ImputeMethodSample,
                             maxMissingRate: float, afMismatchPthresh: float,
                             minGtForInternalImput: int, ignoreFilterField: bool) =
  ## Compute polygenic scores.
  ##
  ##   scores: seq[float] to which the scores will be written. Will be resized
  ##           to the number of samples in genotypeVcf.
  ##   scoreFile: an open ScoreFile describing the polygenic score.
  ##   genotypeVcf: an open VCF containing genotypes of samples for which to
  ##                calculate scores.
  ##   restrictToCoveredRgns: true if analysis should be restricted to
  ##                          well-called regions in coveredIvals; false to
  ##                          consider all loci to be well-covered.
  ##   coveredIvals: object containing genome regions which have been well-
  ##                 called (covered) by the genotyping method.
  ##   imputeMethodLocus: Imputation method to use when a whole locus fails or
  ##                      is missing / not covered.
  ##   imputeMethodSample: Imputation method to use for individual samples with
  ##                       missing genotype, in a locus that passes QC filters.
  ##   maxMissingRate: loci with more than this rate of missing samples fail
  ##                   QC and are imputed.
  ##   afMismatchPthresh: p-value threshold to warn about allele frequency
  ##                      mismatch between the cohort in genotypeVcf and the
  ##                      polygenic score in scoreFile.
  ##   minGtForInternalImput: Minimum number of genotyped samples at a locus for
  ##                          internal imputation to be applied.
  let nsamples = genotypeVcf.n_samples

  # Initialise the scores to zero; offset will be added later
  scores.setLen(nsamples)
  for i in 0..scores.high:
    scores[i] = 0.0

  # Iterate over PS loci.  For each locus, get its (possibly imputed) dosages,
  # and add its score contribution to the accumulating scores.
  var nloci = 0
  var dosages = newSeqUninitialized[float](nsamples)
  for scoreEntry in scoreFile:
    if getImputedDosages(dosages, scoreEntry, genotypeVcf, restrictToCoveredRgns,
                         coveredIvals, imputeMethodLocus, imputeMethodMissing, 
                         imputeMethodSample, maxMissingRate, afMismatchPthresh, 
                         minGtForInternalImput, ignoreFilterField) == true:
      for i in 0..scores.high:
        scores[i] += dosages[i] * scoreEntry.beta
      nloci += 1

  # Average over the total ploidy to match PLINK behaviour
  for i in 0..scores.high:
    scores[i] /= nloci.toFloat*2.0

  # Add the offset
  for i in 0..scores.high:
    scores[i] += scoreFile.offset


proc main() =
  let doc = """
  Compute polygenic scores from a VCF/BCF.

  Usage:
    nimpress [options] <scoredef> <genotypes.vcf>
    nimpress (-h | --help)
    nimpress --version

  Options:
    -h --help          Show this screen.
    --version          Show version.
    --cov=<path>       Path to a BED file supplying genome regions that have been
                       genotyped in the genotypes.vcf file.
    --imp-locus=<m>    Imputation to apply for whole loci which are either not
                       in the sequenced BED regions, or fail (too many samples 
                       with missing genotype, as set by --maxmis, or failing VCF
                       QUAL field if --ignorefilt is not set). Valid values are 
                       ps, homref, fail, ignore [default: ps].
    --imp-missing=<m>  Imputation to apply for loci which are in the sequenced BED
                       regions (and thus should have been genotyped), but are 
                       completely missing from the VCF. Valid values are homref,
                       ignore [default: homref].
    --imp-sample=<m>   Imputation to apply for an individual sample with missing 
                       genotype. Valid values are ps, homref, fail, int_fail, 
                       int_ps [default: int_ps].
    --maxmis=<f>       Maximum fraction of samples with missing genotypes allowed
                       at a locus. Loci containing more than this fraction of 
                       samples missing will be considered bad, and have all 
                       genotypes (even non-missing ones) imputed [default: 0.05].
    --mincs=<n>        Minimum number of genotypes.vcf samples without missing 
                       genotype at a locus for this locus to be eligible for 
                       internal imputation [default: 100].
    --afmisp=<f>       p-value threshold for warning about allele frequency 
                       mismatch between the polygenic score and the supplied 
                       cohort [default: 0.001].
    --ignorefilt       Ignore the VCF FILTER field. If set, all variants in the 
                       VCF will be used regardless of FILTER field contents. If
                       not set, variants with a FILTER field other than "." or
                       "PASS" will always be imputed.

  Imputation methods:
  ps        Impute with dosage based on the polygenic score effect allele 
            frequency.
  homref    Impute to homozygous reference genotype.
  fail      Do not impute, but fail. Failed samples will have a score of "nan"
  ignore    Completely ignore missing loci, as if they were never in the score 
            definition.
  int_ps    Impute with dosage calculated from non-missing samples in the 
            cohort. At least --mincs non-missing samples must be available for 
            this method to be used, else it will fall back to ps.
  int_fail  Impute with dosage calculated from non-missing samples in the 
            cohort. At least --mincs non-missing samples must be available for 
            this method to be used, else it will fall back to fail.
  """

  let args = docopt(doc, version = "nimpress 0.0.1")

  var consoleLog = newConsoleLogger()
  addHandler(consoleLog)

  let
    maxMissingRate = parseFloat($args["--maxmis"])
    afMismatchPthresh = parseFloat($args["--afmisp"])
    minInternalImputeCohortSize = parseInt($args["--mincs"])
    imputeMethodLocus = parseEnum[ImputeMethodLocus]($args["--imp-locus"])
    imputeMethodMissing = parseEnum[ImputeMethodMissing]($args["--imp-missing"])
    imputeMethodSample = parseEnum[ImputeMethodSample]($args["--imp-sample"])

  var
    genotypeVcf: VCF
    scoreFile: ScoreFile
    coveredIvals: GenomeIntervals
    restrictToCoveredRgns: bool
    ignoreFilterField: bool

  if not open(genotypeVcf, $args["<genotypes.vcf>"]):
    log(lvlFatal, "Could not open input VCF file " & $args["<genotypes.vcf>"])
    quit(-1)

  if not open(scoreFile, $args["<scoredef>"]):
    log(lvlFatal, "Could not open polygenic score file " & $args["<scoredef>"])
    quit(-1)

  restrictToCoveredRgns = false
  if args["--cov"]:
    restrictToCoveredRgns = true
    if not loadBedIntervals(coveredIvals, $args["--cov"]):
      log(lvlFatal, "Could not open coverage BED file " & $args["--cov"])

  ignoreFilterField = false
  if args["--ignorefilt"]:
    ignoreFilterField = true

  var scores = newSeqUninitialized[float](0)       # Will be resized as needed
  computePolygenicScores(scores, scoreFile, genotypeVcf, restrictToCoveredRgns,
                         coveredIvals, imputeMethodLocus, imputeMethodMissing,
                         imputeMethodSample, maxMissingRate, afMismatchPthresh,
                         minInternalImputeCohortSize, ignoreFilterField)

  for i in 0..scores.high:
    echo $samples(genotypeVcf)[i] & "\t" & $scores[i]


when isMainModule:
  main()



import math
import sequtils
import unittest

import hts

import nimpress


proc isNaN(x:float): bool =
  result = x.classify == fcNaN


proc checkFloats(x: seq[float], target: seq[float]): bool = 
  if x.len != target.len:
    return false
  for (xi, ti) in zip(x, target):
      if ti.isNaN != xi.isNaN:
        return false
      if ti.isNaN == false and abs(ti - xi) > 1e-4:
        return false
  return true


suite "set1":
  setup:
    var genotypeVcf:VCF
    var scoreFile:ScoreFile
    var coveredBed:GenomeIntervals
    var scores = newSeqUninitialized[float](0)
    discard open(genotypeVcf, "tests/set1.vcf.gz")
    discard open(scoreFile, "tests/set1.score")
    discard loadBedIntervals(coveredBed, "tests/set1.bed")


  test "Locus:ps_MMR1 Sample:fail Coverage:ignored":
    computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
                           ImputeMethodLocus.ps,
                           ImputeMethodMissing.homref,
                           ImputeMethodSample.fail,
                           maxMissingRate = 1.0,
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100,
                           ignoreFilterField = false)
    check(checkFloats(scores, @[NaN, 0.108, NaN, NaN, NaN, NaN]))


  test "Locus:ps_MMR.2 Sample:fail Coverage:ignored":
    computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodMissing.homref,
                           ImputeMethodSample.fail,
                           maxMissingRate = 0.2, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100,
                           ignoreFilterField = false)
    check(checkFloats(scores, @[0.075166667, 0.1085, NaN, NaN, NaN, -0.0165]))


  test "Locus:ps_MMR1 Sample:homref Coverage:ignored":
    computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodMissing.homref,
                           ImputeMethodSample.homref,
                           maxMissingRate = 0.2, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100,
                           ignoreFilterField = false)
    check(checkFloats(scores, @[0.075166667, 0.1085, 0.075166667, 0.141833333, 0.000166667, -0.0165]))


  test "Locus:ps_MMR1 Sample:int3_ps Coverage:ignored":
    computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodMissing.homref,
                           ImputeMethodSample.int_ps,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 3,
                           ignoreFilterField = false)
    check(checkFloats(scores, @[0.075166667, 0.108, 0.070166667, 0.036833333, 0.006833333, -0.0165]))


  test "Locus:ps_MMR1 Sample:int100_ps Coverage:ignored":
    computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodMissing.homref,
                           ImputeMethodSample.int_ps,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100,
                           ignoreFilterField = false)
    check(checkFloats(scores, @[0.075166667, 0.108, 0.074333333, 0.140333333, 0.006833333, -0.0165]))


  test "Locus:ps_MMR1 Sample:int100_fail Coverage:ignored":
    computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodMissing.homref,
                           ImputeMethodSample.int_fail,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100,
                           ignoreFilterField = false)
    check(checkFloats(scores, @[NaN, 0.108, NaN, NaN, NaN, NaN]))


  test "Locus:homref_MMR1 Sample:fail Coverage:ignored":
    computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
                           ImputeMethodLocus.homref,
                           ImputeMethodMissing.homref,
                           ImputeMethodSample.fail,
                           maxMissingRate = 1.0,
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100,
                           ignoreFilterField = false)
    check(checkFloats(scores, @[NaN, 0.098, NaN, NaN, NaN, NaN]))


  test "Locus:homref_MMR.2 Sample:fail Coverage:ignored":
    computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
                           ImputeMethodLocus.homref, 
                           ImputeMethodMissing.homref,
                           ImputeMethodSample.fail,
                           maxMissingRate = 0.2, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100,
                           ignoreFilterField = false)
    check(checkFloats(scores, @[0.064666667, 0.098, NaN, NaN, NaN, -0.027]))


  test "Locus:homref_MMR1 Sample:homref Coverage:ignored":
    computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
                           ImputeMethodLocus.homref, 
                           ImputeMethodMissing.homref,
                           ImputeMethodSample.homref,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100,
                           ignoreFilterField = false)
    check(checkFloats(scores, @[0.064666667, 0.098, 0.064666667, 0.131333333, -0.010333333, -0.027]))


  test "Locus:fail_MMR1 Sample:fail Coverage:ignored":
    computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
                           ImputeMethodLocus.fail,
                           ImputeMethodMissing.homref,
                           ImputeMethodSample.fail,
                           maxMissingRate = 1.0,
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100,
                           ignoreFilterField = false)
    check(checkFloats(scores, @[NaN, NaN, NaN, NaN, NaN, NaN]))


  test "Locus:fail_MMR.2 Sample:fail Coverage:ignored":
    computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
                           ImputeMethodLocus.fail, 
                           ImputeMethodMissing.homref,
                           ImputeMethodSample.fail,
                           maxMissingRate = 0.2, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100,
                           ignoreFilterField = false)
    check(checkFloats(scores, @[NaN, NaN, NaN, NaN, NaN, NaN]))


  test "Locus:ps_MMR1 Sample:ps Coverage:filtered":
    computePolygenicScores(scores, scoreFile, genotypeVcf, true, coveredBed, 
                           ImputeMethodLocus.ps,
                           ImputeMethodMissing.homref,
                           ImputeMethodSample.ps,
                           maxMissingRate = 1.0,
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100,
                           ignoreFilterField = false)
    check(checkFloats(scores, @[0.081, 0.081, 0.081, 0.1545, 0.006, 0.006]))


  test "Versus PLINK 1.90 defaults":
    computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
                           ImputeMethodLocus.ignore,
                           ImputeMethodMissing.ignore,
                           ImputeMethodSample.int_ps,
                           maxMissingRate = 1.0,
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 0,
                           ignoreFilterField = true)
    # Offset of 0.123 is as PLINK does not use external offsets.
    check(checkFloats(scores, @[0.123-0.03, 0.123-0.01, 0.123-0.076, 0.123-0.096, 0.123-0.132, 0.123-0.16]))


  # TODO: Do not fully understand this algorithm in PLINK yet, so setting aside
  # for now.
  # test "Versus PLINK 1.90 no-mean-imputation":
  #   computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
  #                          ImputeMethodLocus.ignore,
  #                          ImputeMethodMissing.ignore,
  #                          ImputeMethodSample.homref,
  #                          maxMissingRate = 1.0,
  #                          afMismatchPthresh = 1.0, 
  #                          minGtForInternalImput = 0,
  #                          ignoreFilterField = true)
  #   check(checkFloats(scores, @[0.123-0.0375, 0.123-0.01, 0.123-0.15, 0.123+0.025, 0.123-0.4, 0.123-0.3]))


  # test "Versus PLINK 2.00":
  #   computePolygenicScores(scores, scoreFile, genotypeVcf, false, coveredBed, 
  #                          ImputeMethodLocus.ps,
  #                          ImputeMethodMissing.ignore,
  #                          ImputeMethodSample.ps,
  #                          maxMissingRate = 1.0,
  #                          afMismatchPthresh = 1.0, 
  #                          minGtForInternalImput = 0,
  #                          ignoreFilterField = true)
  #   check(checkFloats(scores, @[0.123-0.0294, 0.123-0.01, 0.123-0.884, 0.123+0.0208, 0.123-0.1394, 0.123-0.1674]))

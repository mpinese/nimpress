import math
import unittest

import hts

import nimpress


proc isNaN(x:float): bool =
  result = x.classify == fcNaN


proc checkFloats(x: seq[float], target: seq[float]): bool = 
  result = true
  if x.len != target.len:
    result = false
  else:
    for i in 0..target.high:
      if target[i].isNaN != x[i].isNaN:
        result = false
        break
      if target[i].isNaN == false:
        if abs(target[i] - x[i]) > 1e-4:
          result = false
          break
  if result == false:
    echo("checkFloats failed. Target: " & $target & " != value : " & $x)
  else:
    echo("checkFloats PASS")


suite "set1":
  setup:
    var genotypeVcf:VCF
    var scoreFile:ScoreFile
    var coveredBed:File
    var scores = newSeqUninitialized[float](0)
    discard open(genotypeVcf, "tests/set1.vcf.gz")
    discard open(scoreFile, "tests/set1.score")


  test "Locus:ps_MMR1 Sample:fail":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps,
                           ImputeMethodSample.fail,
                           maxMissingRate = 1.0,
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(checkFloats(scores, @[NaN, 0.093, NaN, NaN, NaN, NaN]))


  test "Locus:ps_MMR.2 Sample:fail":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.fail,
                           maxMissingRate = 0.2, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(checkFloats(scores, @[0.0273333333333333, 0.094, NaN, NaN, NaN, -0.156]))


  test "Locus:ps_MMR1 Sample:homref":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.homref,
                           maxMissingRate = 0.2, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(checkFloats(scores, @[0.0273333333333333, 0.094, 0.0273333333333333, 0.160666666666667, -0.122666666666667, -0.156]))


  test "Locus:ps Sample:int3_ps":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.int_ps,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 3)
    check(checkFloats(scores, @[0.0273333333333333, 0.093, 0.0173333333333333, -0.0493333333333333, -0.109333333333333, -0.156]))


  test "Locus:ps Sample:int100_ps":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.int_ps,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(checkFloats(scores, @[0.0273333333333333, 0.093, 0.0256666666666667, 0.157666666666667, -0.109333333333333, -0.156]))


  test "Locus:ps Sample:int100_fail":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.int_fail,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(checkFloats(scores, @[NaN, 0.093, NaN, NaN, NaN, NaN]))


  test "Locus:homref_MMR1 Sample:fail":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.homref,
                           ImputeMethodSample.fail,
                           maxMissingRate = 1.0,
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(checkFloats(scores, @[NaN, 0.073, NaN, NaN, NaN, NaN]))


  test "Locus:homref_MMR.2 Sample:fail":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.homref, 
                           ImputeMethodSample.fail,
                           maxMissingRate = 0.2, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(checkFloats(scores, @[0.00633333333333334, 0.073, NaN, NaN, NaN, -0.177]))


  test "Locus:homref_MMR1 Sample:homref":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.homref, 
                           ImputeMethodSample.homref,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(checkFloats(scores, @[0.00633333333333334, 0.073, 0.00633333333333334, 0.139666666666667, -0.143666666666667, -0.177]))


  test "Locus:fail_MMR1 Sample:fail":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.fail,
                           ImputeMethodSample.fail,
                           maxMissingRate = 1.0,
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(checkFloats(scores, @[NaN, NaN, NaN, NaN, NaN, NaN]))


  test "Locus:fail_MMR.2 Sample:fail":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.fail, 
                           ImputeMethodSample.fail,
                           maxMissingRate = 0.2, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(checkFloats(scores, @[NaN, NaN, NaN, NaN, NaN, NaN]))


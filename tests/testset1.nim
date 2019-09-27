import unittest

include nimpress



suite "set1":
  setup:
    var genotypeVcf:VCF
    var scoreFile:ScoreFile
    var coveredBed:File
    var scores = newSeqUninitialized[float](0)
    discard open(genotypeVcf, "tests/set1.vcf.gz")
    discard open(scoreFile, "tests/set1.score")


  test "Locus:ps Sample:fail":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps,
                           ImputeMethodSample.fail,
                           maxMissingRate = 1.0,
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[NaN, 0.142333333, NaN, NaN, NaN, NaN])

    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.fail,
                           maxMissingRate = 0.2, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[0.076666667, 0.143333333, 0.093333333, NaN, NaN, -0.04])


  test "Locus:ps Sample:homref":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.homref,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[0.075666667, 0.142333333, 0.092333333, 0.242333333, -0.074333333, -0.041])


  test "Locus:ps Sample:int_ps":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.int_ps,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 3)
    check(scores == @[0.076666667, 0.142333333, 0.093333333, 0.033333333, -0.06, -0.04])

    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.int_ps,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[0.076666667, 0.142333333, 0.093333333, 0.240333333, -0.072666667, -0.04])


  test "Locus:ps Sample:int_fail":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.int_fail,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[0.076666667, 0.142333333, 0.093333333, NaN, NaN, -0.04])



  test "Locus:homref Sample:fail":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps,
                           ImputeMethodSample.fail,
                           maxMissingRate = 1.0,
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[NaN, 0.122333333, NaN, NaN, NaN, NaN])

    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.fail,
                           maxMissingRate = 0.2, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[0.055666667, 0.122333333, 0.072333333, NaN, NaN, -0.061])


  test "Locus:homref Sample:homref":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.homref,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[0.055666667, 0.122333333, 0.072333333, 0.222333333, -0.094333333, -0.061])


  test "Locus:homref Sample:int_ps":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.int_ps,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 3)
    check(scores == @[0.056666667, 0.122333333, 0.073333333, 0.013333333, -0.08, -0.06])

    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.int_ps,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[0.056666667, 0.122333333, 0.073333333, 0.220333333, -0.08, -0.06])


  test "Locus:homref Sample:int_fail":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.int_fail,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[NaN, 0.122333333, NaN, NaN, NaN, NaN])



  test "Locus:fail Sample:fail":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps,
                           ImputeMethodSample.fail,
                           maxMissingRate = 1.0,
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[NaN, NaN, NaN, NaN, NaN, NaN])

    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.fail,
                           maxMissingRate = 0.2, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[NaN, NaN, NaN, NaN, NaN, NaN])


  test "Locus:fail Sample:homref":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.homref,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[NaN, NaN, NaN, NaN, NaN, NaN])


  test "Locus:fail Sample:int_ps":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.int_ps,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 3)
    check(scores == @[NaN, NaN, NaN, NaN, NaN, NaN])

    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.int_ps,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[NaN, NaN, NaN, NaN, NaN, NaN])


  test "Locus:fail Sample:int_fail":
    computePolygenicScores(scores, scoreFile, genotypeVcf, coveredBed, 
                           ImputeMethodLocus.ps, 
                           ImputeMethodSample.int_fail,
                           maxMissingRate = 1.0, 
                           afMismatchPthresh = 1.0, 
                           minGtForInternalImput = 100)
    check(scores == @[NaN, NaN, NaN, NaN, NaN, NaN])



# ./nimpress --maxmis=1 --imp-locus=ps --imp-sample=fail tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=0.2 --imp-locus=ps --imp-sample=fail tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=1 --imp-locus=ps --imp-sample=homref tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=1 --mincs=3 --imp-locus=ps --imp-sample=int_ps tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=1 --mincs=100 --imp-locus=ps --imp-sample=int_ps tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=1 --mincs=3 --imp-locus=ps --imp-sample=int_fail tests/set1.score tests/set1.vcf.gz

# ./nimpress --maxmis=1 --imp-locus=homref --imp-sample=fail tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=0.2 --imp-locus=homref --imp-sample=fail tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=1 --imp-locus=homref --imp-sample=homref tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=1 --mincs=3 --imp-locus=homref --imp-sample=int_ps tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=1 --mincs=100 --imp-locus=homref --imp-sample=int_ps tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=1 --mincs=3 --imp-locus=homref --imp-sample=int_fail tests/set1.score tests/set1.vcf.gz

# ./nimpress --maxmis=1 --imp-locus=fail --imp-sample=fail tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=0.2 --imp-locus=fail --imp-sample=fail tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=1 --imp-locus=fail --imp-sample=homref tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=1 --mincs=3 --imp-locus=fail --imp-sample=int_ps tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=1 --mincs=100 --imp-locus=fail --imp-sample=int_ps tests/set1.score tests/set1.vcf.gz
# ./nimpress --maxmis=1 --mincs=3 --imp-locus=fail --imp-sample=int_fail tests/set1.score tests/set1.vcf.gz

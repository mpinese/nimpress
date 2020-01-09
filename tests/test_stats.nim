import math
import unittest

import nimpress

let REL_ERROR_THRESHOLD = 1e-5

proc relerror(val: float, target: float): float = abs((val-target)/target)

proc check_relerror(val: float, target: float): bool = relerror(val, target) < REL_ERROR_THRESHOLD


suite "stats":
  test "dbinom":
    check:
      dbinom(0, 0, 0.0) == 1.0
      dbinom(0, 1, 0.0) == 1.0
      dbinom(1, 1, 0.0) == 0.0

      dbinom(0, 0, 1.0) == 1.0
      dbinom(0, 1, 1.0) == 0.0
      dbinom(1, 1, 1.0) == 1.0

      check_relerror(dbinom( 4, 38, 0.25577630), 1.372476e-02)
      check_relerror(dbinom(11, 14, 0.38016959), 2.078122e-03)
      check_relerror(dbinom(11, 15, 0.03339329), 6.860927e-14)
      check_relerror(dbinom( 5, 42, 0.20535736), 6.290288e-02)
      check_relerror(dbinom(12, 15, 0.46804884), 7.570433e-03)
      check_relerror(dbinom(13, 41, 0.21516811), 4.225388e-02)
      check_relerror(dbinom( 1,  2, 0.29185904), 4.133547e-01)
      check_relerror(dbinom(10, 20, 0.95069797), 9.455226e-09)
      check_relerror(dbinom( 6,  6, 0.75765083), 1.891536e-01)
      check_relerror(dbinom( 5,  7, 0.10843472), 2.502462e-04)
      check_relerror(dbinom(30, 33, 0.27261212), 2.446398e-14)
      check_relerror(dbinom(12, 13, 0.45709316), 5.871179e-04)
      check_relerror(dbinom( 5,  6, 0.64056583), 2.325892e-01)
      check_relerror(dbinom( 7, 23, 0.04018314), 2.151619e-05)
      check_relerror(dbinom( 6,  6, 0.92844315), 6.405188e-01)
      check_relerror(dbinom(26, 32, 0.16796189), 2.155554e-15)
      check_relerror(dbinom(33, 36, 0.49076694), 5.933799e-08)
      check_relerror(dbinom( 4,  7, 0.64614847), 2.703090e-01)
      check_relerror(dbinom( 1, 17, 0.95302045), 9.122048e-21)
      check_relerror(dbinom( 0,  5, 0.38257322), 8.972786e-02)

      check_relerror(dbinom(   0, 1000, 0.1), 1.747871e-46)
      check_relerror(dbinom( 100, 1000, 0.1), 4.201679e-02)
      check_relerror(dbinom( 200, 1000, 0.1), 1.639377e-21)
      check_relerror(dbinom( 500, 1000, 0.5), 2.522502e-02)
      check_relerror(dbinom(4000, 8000, 0.5), 8.920342e-03)


  test "pbinom":
    check:
      check_relerror(pbinom( 4, 38, 0.25577630), 1.958270e-02)
      check_relerror(pbinom(11, 14, 0.38016959), 9.996500e-01)
      check_relerror(pbinom(11, 15, 0.03339329), 1.000000e+00)
      check_relerror(pbinom( 5, 42, 0.20535736), 1.120903e-01)
      check_relerror(pbinom(12, 15, 0.46804884), 9.982583e-01)
      check_relerror(pbinom(13, 41, 0.21516811), 9.571488e-01)
      check_relerror(pbinom( 1,  2, 0.29185904), 9.148183e-01)
      check_relerror(pbinom(10, 20, 0.95069797), 9.918892e-09)
      check_relerror(pbinom( 6,  6, 0.75765083), 1.000000e+00)
      check_relerror(pbinom( 5,  7, 0.10843472), 9.999897e-01)
      check_relerror(pbinom(30, 33, 0.27261212), 1.000000e+00)
      check_relerror(pbinom(12, 13, 0.45709316), 9.999620e-01)
      check_relerror(pbinom( 5,  6, 0.64056583), 9.309152e-01)
      check_relerror(pbinom( 7, 23, 0.04018314), 9.999981e-01)
      check_relerror(pbinom( 6,  6, 0.92844315), 1.000000e+00)
      check_relerror(pbinom(26, 32, 0.16796189), 1.000000e+00)
      check_relerror(pbinom(33, 36, 0.49076694), 1.000000e+00)
      check_relerror(pbinom( 4,  7, 0.64614847), 4.765519e-01)
      check_relerror(pbinom( 1, 17, 0.95302045), 9.148500e-21)
      check_relerror(pbinom( 0,  5, 0.38257322), 8.972786e-02)

      check_relerror(pbinom(   0, 1000, 0.1), 1.747871e-46)
      check_relerror(pbinom( 100, 1000, 0.1), 5.265991e-01)
      check_relerror(pbinom( 200, 1000, 0.1), 1.000000e+00)
      check_relerror(pbinom( 500, 1000, 0.5), 5.126125e-01)
      check_relerror(pbinom(4000, 8000, 0.5), 5.044602e-01)


  test "binom_test":
    check:
      binom_test(0, 10, 0.0) == 1.0
      binom_test(1, 10, 0.0) == 0.0
      binom_test(10, 10, 0.0) == 0.0

      binom_test(0, 10, 1.0) == 0.0
      binom_test(1, 10, 1.0) == 0.0
      binom_test(10, 10, 1.0) == 1.0

      check_relerror(binom_test( 4, 38, 0.25577630), 3.886479e-02)
      check_relerror(binom_test(11, 14, 0.38016959), 3.663595e-03)
      check_relerror(binom_test(11, 15, 0.03339329), 6.940568e-14)
      check_relerror(binom_test( 5, 42, 0.20535736), 1.861900e-01)
      check_relerror(binom_test(12, 15, 0.46804884), 1.669148e-02)
      check_relerror(binom_test(13, 41, 0.21516811), 1.270991e-01)
      check_relerror(binom_test( 1,  2, 0.29185904), 4.985364e-01)
      check_relerror(binom_test(10, 20, 0.95069797), 9.918892e-09)
      check_relerror(binom_test( 6,  6, 0.75765083), 3.466710e-01)
      check_relerror(binom_test( 5,  7, 0.10843472), 2.605677e-04)
      check_relerror(binom_test(30, 33, 0.27261212), 2.537229e-14)
      check_relerror(binom_test(12, 13, 0.45709316), 9.811251e-04)
      check_relerror(binom_test( 5,  6, 0.64056583), 4.296175e-01)
      check_relerror(binom_test( 7, 23, 0.04018314), 2.345121e-05)
      check_relerror(binom_test( 6,  6, 0.92844315), 1.000000e+00)
      check_relerror(binom_test(26, 32, 0.16796189), 2.255836e-15)
      check_relerror(binom_test(33, 36, 0.49076694), 8.212623e-08)
      check_relerror(binom_test( 4,  7, 0.64614847), 7.038423e-01)
      check_relerror(binom_test( 1, 17, 0.95302045), 9.148500e-21)
      check_relerror(binom_test( 0,  5, 0.38257322), 1.640556e-01)

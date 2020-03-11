import math
import unittest

import nimpress


# Deliberate fail case to test CI's ability to detect failures
suite "fails":
  test "fails":
    check:
      1 == 0

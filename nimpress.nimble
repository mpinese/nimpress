# Package

version       = "0.0.1"
author        = "Mark Pinese <mpinese@ccia.org.au>"
description   = "Calculate polygenic scores from VCFs"
license       = "MIT"
srcDir        = "src"
bin           = @["nimpress"]



# Dependencies

requires "nim >= 1.0.0", "hts >= 0.2.21", "docopt >= 0.6.8"

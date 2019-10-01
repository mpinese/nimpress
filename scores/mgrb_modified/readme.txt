Scores were converted from those used in the Medical Genome Reference Bank analysis,
see https://www.biorxiv.org/content/10.1101/473348v1.  In summary, for each score:
  1. Variants were sourced from the given publications
  2. Only strand-specific biallenic SNPs were retained
  3. SNPs not well-covered by the MGRB sequencing were replaced with high-correlation
     (r2 > 0.9) covered tag SNPs, where this was possible.  If it was not possible,
     uncovered SNPs were dropped.
  4. Scores were transformed to give the risk for the alternate allele only.

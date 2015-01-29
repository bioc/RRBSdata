\name{rrbs}
\alias{rrbs}
\docType{data}
\title{
An RRBS data set with 12 samples and 10,000 simulated DMRs.
}
\description{
To obtain a dataset with known DMRs, we used 12 human control samples
from a previously published RRBS data set
(Schoofs et al., Blood, 2013). We simulated DMRs of different lengths and
intensities by altering the number of methylated reads of some of the
CpG sites in half of the samples, which could be considered as the cancer
samples. For more information on the simulation, see Details.

}
\usage{data(rrbs)}
\format{
  A \code{BSraw-class} object.
}
\details{
We downloaded CpG island positions from UCSC database and
filtered out all islands that were not covered in any of the samples
(27,718 remaining islands). Only islands with not less than 10
covered CpG sites were considered to receive a DMR (24,698 islands). Within these 24,698 CpG
islands we incorporated 10,000 DMRs with methylation differences of
10, 20, 30, or 40\%. The DMRs spanned 10, 20, 30, or
40\% of the CpG sites of the respective islands. We made sure that
we gained the same amount of DMRs for each of
the 16 combinations of methylation difference and percentage of
modified CpG sites, that is, we gained 1,250 DMRs per combination.

We incorporated the DMRs as follows:
1.) We devided the 12 control samples into 6 ``cancer'' samples
and 6 ``control'' samples.
2.) Within each CpG island we only considered regions to receive a
DMR if each of its CpG sites was covered in at least half of all
samples. Those regions are referred to as ``covered island regions''.
3.) For each CpG site within covered island regions we determined its
minimum and maximum smoothed methylation level across all samples.
4.) For each CpG island we determined the maximum percentage of
neighbored CpG sites that can be altered by determining the percentage of CpG sites within
its biggest covered island region on all CpG sites
within the island.
5.) Each of the 10,000 DMRs was sampled into a covered island region
of a CpG island that was appropriate to harbor the DMR, in terms of minimum and maximum
DNA methylation and CpG percentage of its largest covered island
region on all CpG sites within the island. Not more than one DMR was incorporated per CpG island.
6.) Whenever it was possible to increase the DNA methylation by the amount of
the difference of the respective DMR (that is, the resulting
methylation level is below 1), it was increased in the cancer
samples. Otherwise, the DNA methylation was decreased by the amount of
difference.
7.) To increase the DNA methylation in a cancer sample within a region by a certain
amount, the number of methylated reads of the respective CpG sites was
increased. For instance, a CpG site of coverage 12 with 3 methylated
reads (and a relative methylation level of 0.25) within a region
that should get a methylation difference of 0.3, received 4 additional
methylated reads (because: 0.3x12 = 4).

The scripts to simulate the data can be found under /R.
}
\source{
Schoofs T, Rohde C, Hebestreit K, Klein HU, Goellner S, Schulze I,
  Lerdrup M, Dietrich N, Agrawal-Singh S, Witten A, Stoll M, Lengfelder
  E, Hofmann WK, Schlenke P, Buechner T, Hansen K, Berdel WE, Rosenbauer
  F, Dugas M, Mueller-Tidow C (2012). DNA methylation changes are a late
  event in acute promyelocytic leukemia and coincide with loss of
  transcription factor binding. Blood.
}
\references{
Hebestreit K, Klein HU. Global test and BiSeq are the methods of choice for testing genomic regions for differential methylation in bisulfite sequencing data. In preparation.
}
\examples{
data(rrbs)
rrbs
}
\keyword{datasets}

\name{islands}
\alias{islands}
\docType{data}
\title{
CpG islands with information on DMRs.
}
\description{
The 27,718 CpG islands were dowloaded from UCSC database. Only islands with not less than 10
covered CpG sites were considered to receive a DMR (24,698
islands). Within these islands we incorporated 10,000 DMRs.
}
\usage{data(islands)}
\format{
  A \code{GRanges-class} object with metdata columns:

  \itemize{
    \item \code{cluster.id} An ID of the CpG island.
    \item \code{cpg.no} Number of CpG sites within CpG island.
    \item \code{cpg.perc} Percentage of CpG sites on all nucleotides
    within island
    (\code{cpg.no/width(islands)}).
    \item \code{ave.cov} Average coverage of CpG sites within CpG island.
    \item \code{KStest.clust} Test statistic of test on complete spatial randomness of CpG
    sites within CpG island.
    \item \code{KStest.pval.clust} P value of test on complete spatial
    randomness of CpG sites within CpG island.
    \item \code{dmr.start} Start site position of the incorporated
    DMR. Is \code{NA}, if no DMR is incorporated for this CpG island.
    \item \code{dmr.end} End site position of the incorporated
      DMR. Is \code{NA}, if no DMR is incorporated for this CpG island.
    \item \code{dmr.meth.diff} Methylation difference of the
    incorporated DMR. Is \code{NA}, if no DMR is incorporated for this
    CpG island.
    \item \code{cpg.no.dmr} Number of CpG sites within incorporated
    DMR. Is \code{NA}, if no DMR is incorporated for this CpG island.
    \item \code{cpg.perc.dmr} Percentage of CpG sites within
    incorporated DMR on all CpG sites within this island. Is \code{NA},
    if no DMR is incorporated for this CpG island.
    \item \code{ave.cov.dmr} The average coverage of all CpG sites
    within incorporated DMR. Is \code{NA}, if no DMR is incorporated for
    this CpG island.
    \item \code{KStest.clust.dmr} Test statistic of test on complete spatial randomness of CpG
    sites within incorporated DMR. Is \code{NA}, if no DMR is
    incorporated for this CpG island.
    \item \code{KStest.pval.clust.dmr} P value of test on complete spatial randomness of CpG
    sites within incorporated DMR. Is \code{NA}, if no DMR is
    incorporated for this CpG island.
  }
}
\references{
Hebestreit K, Klein HU. Global test and BiSeq are the methods of choice for testing genomic regions for differential methylation in bisulfite sequencing data. In preparation.
}
\examples{
data(islands)
head(islands)
}
\keyword{datasets}

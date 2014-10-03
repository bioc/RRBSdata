\name{diffMethCpGs}
\alias{diffMethCpGs}
\docType{data}
\title{
The differentially methylated CpG sites.
}
\description{
This object comprises all CpG sites that are differentially methylated
between 'cancer' and 'control' samples.
}
\usage{data(diff.meth.cpgs)}
\format{ A \code{GRanges-class} object with information on the CpG
  island (\code{cluster.id}), the covered island region
  (\code{cluster.id.2}) and the simulated methylation difference (\code{meth.diff}).
}
\references{
Hebestreit K, Klein HU. Global test and BiSeq are the methods of choice for testing genomic regions for differential methylation in bisulfite sequencing data. In preparation.
}
\examples{
data(diffMethCpGs)
head(diffMethCpGs)
}
\keyword{datasets}

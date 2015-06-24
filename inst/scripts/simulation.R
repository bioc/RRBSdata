### Author: Katja Hebestreit (katjah@stanford.edu)
### Date: January 29th 2015

library(BiSeq)
library(spatstat)
source("R/simulateDMRs.R")

# What you need:
# rrbs: a BSraw object
# islands: a GRanges with CpG islands or any other genomic regions where you want to have DMRs
# dmr.attr: a data frame with the attributes of the DMRs you want to simulate; rows correspond to DMRs

# Please check the simulateDMRs.R file for all input and output objects.

## Define DMR attributes:
abs.meth.diff <- c(0.1, 0.2, 0.3, 0.4)
CpG.perc <- c(0.1, 0.2, 0.3, 0.4)
n <- 2500
dmr.attr <- data.frame(CpG.perc = rep(abs.meth.diff, n),
                       abs.meth.diff = rep(CpG.perc, each = n))

# optional:
set.seed(1234)

## Simulate DMRs (please see the simulateDMRs.R file to understand the output):
sim <- simulateDMRs(object = rrbs,
                    ind.samples = which(colData(rrbs)$group == "cancer"),
                    regions = islands,
                    dmr.attr = dmr.attr,
                    h=80,
                    mc.cores=12,
                    perc.samples = 0.5,
                    min.no.cpg = 10)

# takes approx. 20h
# quite often:
# Warning message:
# In ks.test(U, "punif", ...) :
# ties should not be present for the Kolmogorov-Smirnov test
# --> That does not effect the simulation.


## This is all optional:

rrbs.dmrs <- sim$rrbs.dmrs
regions.dmrs <- sim$regions[!is.na(sim$regions$cpg.perc)]
regions.dmrs <- sort(regions.dmrs)

## A possible test that it worked:
plotMethMap(rrbs.dmrs, regions.dmrs[2], groups = colData(rrbs.dmrs)$group)


## Add region information to regions without DMRs:

rrbs <- sim$rrbs
regions <- sim$regions
ind.no.dmr <- which(is.na(regions$cpg.no))

for(i in ind.no.dmr){
    ov.rrbs <- rrbs[which(rowRanges(rrbs)$cluster.id == regions$cluster.id[i]),]
    ov.cpg <- rowRanges(ov.rrbs)
    cpg.no <- length(ov.cpg)
    regions[i]$cpg.no <- cpg.no
    regions[i]$cpg.perc <- cpg.no / width(regions[i])
    regions[i]$ave.cov <- median(totalReads(ov.rrbs))
}

sim$regions <- regions


# number of diff. CpG sites per region:
t <- table(sim$diff.meth.cpgs$cluster.id)
df <- as.data.frame(t)
quantile(df)
#  0%  25%  50%  75% 100%
#   1    7   14   30  298

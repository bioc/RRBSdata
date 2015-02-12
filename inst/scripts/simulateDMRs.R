### Author: Katja Hebestreit (katjah@stanford.edu)
### Date: January 29th 2015

# Incorporates DMRs within a given BSraw object

# Input:
# object: BSraw object, colData: group (two groups)
# ind.samples: vector, that specifies the samples, that should get differential methylation
# regions: GRanges of regions within which DMRs should lie
# dmr.attr: data.frame. per DMR one row. columns: CpG.perc (of respective region), abs.meth.diff
# min.no.cpg: minimum number of covered CpG sites within regions
# h: band width for smoothing to determine original min and max of methylation within regions, see predictMeth
# mc.cores: passed to predictMeth
# perc.samples: each CpG site within a DMR should be smoothed in perc.samples*100 per cent of the samples -> covered regions

# Output: a list of
# - rrbs.dmrs: a BSraw object with incorporated DMRs
# - diff.meth.cpgs: a GRanges of diff. meth. cpg sites, rownames= region.id, elementMetadata: region.id, meth.diff
# - regions: input regions (e.g. promoters), rownames: region.id,
#   elementMetadata:
#   region.id, cpg.no, cpg.perc, ave.cov (mean coverage of covered CpG sites),
#   KStest.clust, KStest.pval.clust (test on complete spatial randomness of CpG sites within region)
#   dmr.start, dmr.end, dmr.meth.diff, cpg.no.dmr, cpg.perc.dmr, ave.cov.dmr (mean coverage of covered CpG sites within DMR),
#   KStest.clust.dmr, KStest.pval.clust.dmr (test on complete spatial randomness of DMR-CpG-sites within region)

# Simulates per region not more than one DMR.
# The input regions should not overlap.


simulateDMRs <- function(object, ind.samples, regions, dmr.attr, h=80, mc.cores=6, perc.samples = 0.5, min.no.cpg){
  strand(regions) <- "*"
  if( length(regions) != length(reduce(regions)) ){
    stop("Your regions are not disjunct. But they should be :-)\n")
  }

  strand(object) <- "*"
  object <- sort(object)
  regions <- sort(regions)
  methReads.all.new <- methReads(object)

  diff.meth.cpgs <- GRanges()
  diff.meth.cpgs$cluster.id <- character()
  diff.meth.cpgs$meth.diff <- numeric()

  # cluster.id: equivalent to region id

  regions$cluster.id <- paste("region_", seq(along=regions), sep="")
  names(regions) <- regions$cluster.id
  regions$cpg.no <- as.integer(rep(NA, length=length(regions)))
  regions$cpg.perc <- as.numeric(rep(NA, length=length(regions)))
  regions$ave.cov <- as.numeric(rep(NA, length=length(regions)))
  regions$KStest.clust <- as.numeric(rep(NA, length=length(regions)))
  regions$KStest.pval.clust <- as.numeric(rep(NA, length=length(regions)))
  regions$dmr.start <- as.integer(rep(NA, length=length(regions)))
  regions$dmr.end <- as.integer(rep(NA, length=length(regions)))
  regions$dmr.meth.diff <- as.numeric(rep(NA, length=length(regions)))
  regions$cpg.no.dmr <- as.integer(rep(NA, length=length(regions)))
  regions$cpg.perc.dmr <- as.numeric(rep(NA, length=length(regions)))
  regions$ave.cov.dmr <- as.numeric(rep(NA, length=length(regions)))
  regions$KStest.clust.dmr <- as.numeric(rep(NA, length=length(regions)))
  regions$KStest.pval.clust.dmr <- as.numeric(rep(NA, length=length(regions)))

  overlap <- findOverlaps(rowData(object), regions)

  rowData(object)$cluster.id <- rep(NA, nrow(object))
  rowData(object)$cluster.id[overlap@queryHits] <- regions$cluster.id[overlap@subjectHits]

  object.regions <- object[overlap@queryHits,]
  elementMetadata(rowData(object.regions))$cluster.id <-
    elementMetadata(regions)$cluster.id[overlap@subjectHits]


  cat("Begin to smooth the data within all regions to obtain minimum and maximum relative methylation values.\n")
  object.smoothed <- predictMeth(object.regions, h = h, mc.cores = mc.cores)

  n.s <- ncol(object.smoothed)
  rowData(object.smoothed)$perc.samples <-
    apply(methLevel(object.smoothed), 1, function(x){
      sum(!is.na(x)) / n.s
    })

  # ind.covered: per CpG: is CpG covered in >= perc.samples of samples?
  ind.covered <- rowData(object.smoothed)$perc.samples >= perc.samples

  # cluster.id.2: regions of CpG sites covered in >= perc.samples in respective region (cluster.id)
  cluster.id.2 <- character(length=nrow(object.smoothed))
  cluster.id.2[!ind.covered] <- NA

  # assign cluster.id.2 per CpG site covered in >= perc.samples:
  for( clus in unique(rowData(object.smoothed)$cluster.id)){
    ind.clus <- rowData(object.smoothed)$cluster.id == clus
    if( any(!is.na(cluster.id.2[ind.clus])) ){
      if( all(!is.na(cluster.id.2[ind.clus])) ){
        cluster.id.2[ind.clus] <- paste(clus, "_", 1, sep="")
      } else {
        ind.na <- which(is.na(cluster.id.2[ind.clus]))
        s <- integer(length=sum(ind.clus))
        s[ind.na] <- 1
        s <- cumsum(s)
        s <- s + 1
        ind.1 <- s[-ind.na]
        ind.2 <- length(unique(ind.1))
        ind.3 <- rep(1:ind.2, table(ind.1))
        cluster.id.2[ind.clus][-ind.na] <- paste(clus, "_", ind.3, sep="")
      }
    }
  }
  rowData(object.smoothed)$cluster.id.2 <- cluster.id.2

  ind.ov <- findOverlaps(rowData(object), rowData(object.smoothed))
  rowData(object)$cluster.id.2 <- rep(NA, nrow(object))
  rowData(object)$cluster.id.2[ind.ov@queryHits] <- rowData(object.smoothed)$cluster.id.2[ind.ov@subjectHits]

  # extremes.cpg: data frame with minimum and maximum smoothed methylation value per CpG site
  extremes.cpg <- data.frame(min = apply(methLevel(object.smoothed), 1, min, na.rm=TRUE),
                             max = apply(methLevel(object.smoothed), 1, max, na.rm=TRUE))

  # cov.region.attr: for each cluster of covered CpG sites by >= perc.samples within regions (cluster.id.2):
  # rownames: cluster.id.2
  # - minimum (min.meth) and maximum (max.meth) methylation value
  # - number of covered CpG sites (no.cov.cpg) by >= perc.samples
  # - cluster.id a.k.a. region
  # - number of covered CpG sites (by at least one sample) in respective region (no.all.cpg)
  cov.region.attr <- data.frame(
                       min.meth = tapply(extremes.cpg$max,
                         rowData(object.smoothed)$cluster.id.2,
                         min),
                       max.meth = tapply(extremes.cpg$max,
                         rowData(object.smoothed)$cluster.id.2,
                         max),
                       no.cov.cpg = as.integer(table(rowData(object.smoothed)$cluster.id.2)))
  rownames(cov.region.attr) <- levels(as.factor(rowData(object.smoothed)$cluster.id.2))
  cluster.id <- character(length=nrow(cov.region.attr))
  no.all.cpg <- integer(length=nrow(cov.region.attr))
  for(i in seq(along=rownames(cov.region.attr))) {
    cov.reg <- rownames(cov.region.attr)[i]
    ind.cov.reg <- which(rowData(object.smoothed)$cluster.id.2 == cov.reg)
    cluster.id[i] <- rowData(object.smoothed)$cluster.id[ ind.cov.reg[1] ]
    ind.reg <- which(rowData(object.smoothed)$cluster.id == cluster.id[i])
    no.all.cpg[i] <- length(ind.reg)
  }
  cov.region.attr$cluster.id <- cluster.id
  cov.region.attr$no.all.cpg <- no.all.cpg

  cov.region.attr <- cov.region.attr[cov.region.attr$no.all.cpg >= min.no.cpg,]

  cluster.ids <- unique(cov.region.attr$cluster.id)
  # region.attr: data frame of regions (cluster.id) corresponding to covered CpG clusters by >= perc.samples
  region.attr <- data.frame(min.meth = numeric(length=length(cluster.ids)),
                            max.meth = numeric(length=length(cluster.ids)),
                            max.cov.cpg.perc = numeric(length=length(cluster.ids)))
  # max.cov.cpg.perc : maximum percentage of CpG sites within a cluster within this region
  rownames(region.attr) <- cluster.ids
  for(clus in cluster.ids) {
    part <- cov.region.attr[cov.region.attr$cluster.id == clus,]
    region.attr[clus,"min.meth"] <- min(part$min.meth)
    region.attr[clus,"max.meth"] <- max(part$max.meth)
    region.attr[clus,"max.cov.cpg.perc"] <- max(part$no.cov.cpg) / part$no.all.cpg[1]
  }

  dmr.attr <- dmr.attr[order(dmr.attr$CpG.perc, dmr.attr$abs.meth.diff, decreasing=TRUE),]
  # dmr.attr.table: summarization of DMR characteristics, how many DMRs per characteristic combination
  dmr.attr.table <- unique(dmr.attr)
  dmr.attr.table$no.dmr <- apply(dmr.attr.table,
                                 1,
                                 function(y) {
                                   sum(apply(dmr.attr, 1, identical, y))
                                 }
                                 )

  # per DMR: sample region for this DMR
  dmr.region.assign <- dmr.attr
  dmr.region.assign$region <- character(length=nrow(dmr.attr))

  region.attr.i <- region.attr
  cat("Sampling DMR - region assignment...")
  for(i in 1:nrow(dmr.attr.table)){
    dmr.i <- dmr.attr.table[i,]
    ind.meth.diff <- which(( region.attr.i$min.meth - dmr.i$abs.meth.diff ) >= 0  | ( region.attr.i$max.meth + dmr.i$abs.meth.diff ) <= 1  )
    ind.perc.cov <- which(( region.attr.i$max.cov.cpg.perc ) >= dmr.i$CpG.perc)
    ind.possible <- intersect(ind.meth.diff, ind.perc.cov)
    if(length(ind.possible) < dmr.i$no.dmr) {
      stop("There are no regions left that can take ", dmr.i$no.dmr, " DMR(s) with an absolute methylation difference of ", dmr.i$abs.meth.diff )
    }
    s <- sample(x = rownames(region.attr.i)[ind.possible],
                size = dmr.i$no.dmr)
    dmr.region.assign$region[(dmr.region.assign$CpG.perc == dmr.i$CpG.perc)
                             &
                             (dmr.region.assign$abs.meth.diff == dmr.i$abs.meth.diff)] <- s
    region.attr.i <- region.attr.i[! is.element(rownames(region.attr.i), s),]
  }
  cat(" Done.\n")

  # incorporate DMRs in clusters (cluster.id.2) within regions (cluster.id):
  # per region one DMR only
  cat("Begin to incorporate DMRs... ")
  n.dmrs <- nrow(dmr.region.assign)
  i.dmr <- 0
  pb <- txtProgressBar(min=0, max=n.dmrs, style=3)

  for(i in 1:n.dmrs) {
    i.dmr <- i.dmr + 1
    setTxtProgressBar(pb, value = i.dmr)

    cpg.perc <- dmr.region.assign$CpG.perc[i]
    meth.diff <- dmr.region.assign$abs.meth.diff[i]
    region <- dmr.region.assign$region[i]
    # possible cluster.id.2 within chosen regions:
    possible.clusters <- cov.region.attr[cov.region.attr$cluster.id == region,]
    ind.meth.diff <- which(( possible.clusters$min.meth - meth.diff ) >= 0  | ( possible.clusters$max.meth + meth.diff ) <= 1  )
    ind.perc.cov <- which(( possible.clusters$no.cov.cpg / possible.clusters$no.all.cpg ) >= cpg.perc)
    ind.possible <- intersect(ind.meth.diff, ind.perc.cov)
    if(length(ind.possible) > 1) {
      s <- sample(ind.possible, size=1)
    }
    if(length(ind.possible) == 1){
      s <- ind.possible
    }
    cluster <- possible.clusters[s,]

    ind.region <- which(is.element(rowData(object)$cluster.id, cluster$cluster.id))
    object.region <- object[ind.region,]
    no.cpgs.region <- cluster$no.all.cpg
    no.cpgs <- round(cpg.perc * no.cpgs.region)
    no.cpgs.cov.region <- cluster$no.cov.cpg
    start.cpgs <- sample(1 : (no.cpgs.cov.region - no.cpgs + 1), 1)
    end.cpgs <- start.cpgs + no.cpgs - 1
    ind.cpgs <- which(is.element(rowData(object.region)$cluster.id.2, rownames(cluster)))[start.cpgs : end.cpgs]
    object.dmr <- object.region[ind.cpgs,]
    totalReads <- totalReads(object.dmr)[, ind.samples]
    methReads <- methReads(object.dmr)[, ind.samples]

    summable <- cluster$max.meth + meth.diff <= 1
    if(!summable) {
      meth.diff <- - meth.diff
    }

    add <- meth.diff * totalReads
    dig <- strsplit(as.character(format(add, digits=7)), ".", fixed=TRUE)
    prob <- as.numeric(paste("0.",sapply(dig, function(x) x[2]), sep=""))
    prob[is.na(prob)] <- 0
    sample.add <- as.numeric(sapply(dig, function(x) x[1]))

    if(summable){
      sample.add <- data.frame(add1 = sample.add, add2 = sample.add+1 , prob=prob)
    } else {
      sample.add <- data.frame(add1 = sample.add, add2 = sample.add-1 , prob=prob)
    }
    s.add <- apply(sample.add, 1, function(x){
      sample(x[c(1,2)], size=1, prob=c(1-x[3], x[3]))
    })

    methReads.new <- methReads + s.add

    if(summable) {
      ind.too.much <- methReads.new > totalReads
      methReads.new[ind.too.much] <- totalReads[ind.too.much]
    } else {
      ind.neg <- methReads.new < 0
      methReads.new[ind.neg] <- 0
    }

    storage.mode(methReads.new) <- "integer"

    methReads.all.new[ind.region[ind.cpgs], ind.samples] <- methReads.new

    rowData(object.dmr)$meth.diff <- meth.diff
    diff.meth.cpgs <- c(diff.meth.cpgs, rowData(object.dmr))

    ## test on complete spatial randomness of CpG sites within region / DMR:

    pos.region <- start(object.region)
    start.region <- min(pos.region)
    end.region <- max(pos.region)
    window.region <- owin(xrange=c(start.region, end.region), yrange=c(c(0, 2)))
    ppp.region <- ppp(x = pos.region,
                      y = rep(1, length(pos.region)),
                      window = window.region)
    kstest.region <- kstest(ppp.region,
                            covariate = "x")

    pos.dmr <- start(object.dmr)
    ppp.dmr <- ppp(x = pos.dmr,
                   y = rep(1, length(pos.dmr)),
                   window = window.region)
    kstest.dmr <- kstest(ppp.dmr,
                         covariate = "x")

    # quite often (without duplicate positions):
    # Warning message:
    # In ks.test(U, "punif", ...) :
    # ties should not be present for the Kolmogorov-Smirnov test


    regions[cluster$cluster.id]$cpg.no <- no.cpgs.region
    regions[cluster$cluster.id]$cpg.perc <- no.cpgs.region / width(regions[cluster$cluster.id])
    regions[cluster$cluster.id]$ave.cov <- median(totalReads(object.region))
    regions[cluster$cluster.id]$KStest.clust <- kstest.region$statistic
    regions[cluster$cluster.id]$KStest.pval.clust <- kstest.region$p.value
    regions[cluster$cluster.id]$dmr.start <- start(object.dmr)[1]
    regions[cluster$cluster.id]$dmr.end <- rev(start(object.dmr))[1]
    regions[cluster$cluster.id]$dmr.meth.diff <- meth.diff
    regions[cluster$cluster.id]$cpg.no.dmr <- no.cpgs
    regions[cluster$cluster.id]$cpg.perc.dmr <- cpg.perc
    regions[cluster$cluster.id]$ave.cov.dmr <- median(totalReads(object.dmr))
    regions[cluster$cluster.id]$KStest.clust.dmr <- kstest.dmr$statistic
    regions[cluster$cluster.id]$KStest.pval.clust.dmr <- kstest.dmr$p.value
  }
  close(pb)

  rrbs.dmrs <- object
  methReads(rrbs.dmrs) <- methReads.all.new

  names(diff.meth.cpgs) <- diff.meth.cpgs$cluster.id

  return(
    list(rrbs.dmrs = rrbs.dmrs,
         diff.meth.cpgs = diff.meth.cpgs,
         regions = regions
        )
    )

}

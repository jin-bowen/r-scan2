# Spatial sensitivity applies only to VAF-based calling! Not to mutation signature-based rescue!
# Use small tiles for parallelization (~1MB) because of basepair resolution.
DEPRECATED_compute.spatial.sensitivity.depth <- function(single.cell.id, bulk.id,
    static.filter.params, joint.dptab.path, genome.string, sens.tilewidth=1e3,
    grs.for.sens=genome.string.to.tiling(genome.string, tilewidth=sens.tilewidth, group='auto'),
    grs.for.parallelization=genome.string.to.tiling(genome.string, tilewidth=1e6, group='auto'),
    quiet=TRUE, report.mem=TRUE)
{
    warning("THIS METHOD IS DEPRECATED! IT IS ~10x SLOWER THAN AN EQUIVALENT METHOD USING bedtools map!")
    cat('Gathering read depth data for spatial somatic calling sensitivity using', length(grs.for.parallelization), 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

    sfp <- static.filter.params

    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs.for.parallelization))
        xs <- future.apply::future_lapply(1:length(grs.for.parallelization), function(i) {
            gr <- grs.for.parallelization[i,]

            # Only retain grs.for.sens windows that *start* in the
            # grs.for.parallelization window. This avoids assigning a grs.for.sens
            # window to 2 grs.for.parallelization windows if it spans the boundary.
            # There's probably a better way to do this.
            grs.for.sens2 <- GRanges(seqnames=seqnames(grs.for.sens),
                ranges=IRanges(start=start(grs.for.sens), width=1),
                seqinfo=seqinfo(grs.for.sens))
            gr2 <- grs.for.sens[IRanges::countOverlaps(grs.for.sens2, gr, minoverlap=1) > 0,]

            pc <- perfcheck(paste('read.depth.2sample', i),
                dp <- read.depth.2sample(path=joint.dptab.path, sc.sample=single.cell.id,
                    bulk.sample=bulk.id, keep.coords=TRUE,
                    # N.B. read.depth.2sample requires region to be a GRanges with only
                    # one range in it. Cases arise in regular use where excluding ranges
                    # from grs.for.sens that do not *start* in gr can lead to disjoint
                    # ranges. To be safe, read in depth info for the union of gr2 and gr
                    # (because gr2 can also extend outside of gr). This makes sure all
                    # possible depth data is available AND that only one contigous range
                    # is present in `region`.
                    region=reduce(c(gr, gr2)), quiet=quiet),
                report.mem=report.mem)
            p(class='sticky', pc, amount=0)

            # Set up several GRanges objects for the binnedAverage() method
            pc <- perfcheck(paste("compute.averages", i), {
                    dp.gr <- GRanges(seqnames=dp$chr, ranges=IRanges(dp$pos, width=1),
                        sc.dp=dp[[3]],
                        bulk.dp=dp[[4]],
                        # these are logical vectors: 1 if each single base is >=
                        # the filter cutoff. The binned average can then be multiplied by the width
                        # to recover number of bases.
                        base.gt.snv.sc.min.dp=dp[[3]] >= sfp$snv$min.sc.dp,
                        base.gt.snv.bulk.min.dp=dp[[4]] >= sfp$snv$min.bulk.dp,
                        base.gt.indel.sc.min.dp=dp[[3]] >= sfp$indel$min.sc.dp,
                        base.gt.indel.bulk.min.dp=dp[[4]] >= sfp$indel$min.bulk.dp,
                        seqinfo=seqinfo(gr2))
                    sc.dp <- mcolAsRleList(dp.gr, 'sc.dp')
                    bulk.dp <- mcolAsRleList(dp.gr, 'bulk.dp')
                    base.gt.snv.sc.min.dp <- mcolAsRleList(dp.gr, 'base.gt.snv.sc.min.dp')
                    base.gt.snv.bulk.min.dp <- mcolAsRleList(dp.gr, 'base.gt.snv.bulk.min.dp')
                    base.gt.indel.sc.min.dp <- mcolAsRleList(dp.gr, 'base.gt.indel.sc.min.dp')
                    base.gt.indel.bulk.min.dp <- mcolAsRleList(dp.gr, 'base.gt.indel.bulk.min.dp')
                    gr2 <- binnedAverage(gr2, sc.dp, 'mean.sc.dp')
                    gr2 <- binnedAverage(gr2, bulk.dp, 'mean.bulk.dp')
                    gr2 <- binnedAverage(gr2, base.gt.snv.sc.min.dp, 'mean.gt.snv.sc.min.dp')
                    gr2 <- binnedAverage(gr2, base.gt.snv.bulk.min.dp, 'mean.gt.snv.bulk.min.dp')
                    gr2 <- binnedAverage(gr2, base.gt.indel.sc.min.dp, 'mean.gt.indel.sc.min.dp')
                    gr2 <- binnedAverage(gr2, base.gt.indel.bulk.min.dp, 'mean.gt.indel.bulk.min.dp')
                }, report.mem=report.mem)
            p(class='sticky', pc, amount=1)

            data.table(chr=as.character(seqnames(gr2)), start=start(gr2), end=end(gr2),
                mean.sc.dp=gr2$mean.sc.dp, mean.bulk.dp=gr2$mean.bulk.dp,
                bases.gt.snv.sc.min.dp=width(gr2) * gr2$mean.gt.snv.sc.min.dp,
                bases.gt.snv.bulk.min.dp=width(gr2) * gr2$mean.gt.snv.bulk.min.dp,
                bases.gt.indel.sc.min.dp=width(gr2) * gr2$mean.gt.indel.sc.min.dp,
                bases.gt.indel.bulk.min.dp=width(gr2) * gr2$mean.gt.indel.bulk.min.dp
            )
        })
    }, enable=TRUE)

    return(rbindlist(xs))  # decide if we want to keep this as a GRanges object or not
}



# Spatial sensitivity applies only to VAF-based calling! Not to mutation signature-based rescue!
# IMPORTANT: grs.for.sens must match the GRanges used for compute.spatial.sensitivity.depth().
# N.B. do NOT try to rewrite this function to use a SCAN2 object as argument. future()'s multi-
# core implementation copies the entire parent memory space to each child thread no matter how I
# try to prevent it. For human genomes, this means a waste of ~2-2.5G of RAM per thread.
compute.spatial.sensitivity.abmodel <- function(
    single.cell.id, ab.fits, integrated.table.path, genome.string, sens.tilewidth=1e3,
    grs.for.sens,  # no longer suggest a default. these bins are externally generated
    #grs.for.sens=genome.string.to.tiling(genome.string, tilewidth=sens.tilewidth, group='auto'),
    grs.for.parallelization=analysis.set.tiling.for.parallelization(regions=GenomicRanges::reduce(grs.for.sens), total.tiles=10),
    quiet=TRUE, report.mem=TRUE)
{
    cat('Estimating genome-wide AB for spatial somatic calling sensitivity using', length(grs.for.parallelization), 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')
    # Necessary for AB estimation (infer.gp uses matrix multiplication)
    cat('Setting OpenBLAS corecount to 1. This prevents multithreaded matrix multiplication in chunks where it is undesired.\n')
    RhpcBLASctl::blas_set_num_threads(1)

    cat('Estimating AB at window mid-points..\n')
    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs.for.parallelization))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        xs <- future.apply::future_lapply(1:length(grs.for.parallelization), function(i) {
            gr <- grs.for.parallelization[i,]

            # Only retain grs.for.sens windows that *start* in the
            # grs.for.parallelization window. This avoids assigning a grs.for.sens
            # window to 2 grs.for.parallelization windows if it spans the boundary.
            # There's probably a better way to do this.
            grs.for.sens2 <- GRanges(seqnames=seqnames(grs.for.sens),
                ranges=IRanges(start=start(grs.for.sens), width=1),
                seqinfo=seqinfo(grs.for.sens))
            gr2 <- grs.for.sens[IRanges::countOverlaps(grs.for.sens2, gr, minoverlap=1) > 0,]

            pc <- perfcheck(paste('get.training.sites',i),
                    training.sites <- get.training.sites.for.abmodel.by.range(
                        # N.B. get.training.sites.for.abmodel.by.rnage requires `region` to be
                        # a GRanges with only
                        # one range in it. Cases arise in regular use where excluding ranges
                        # from grs.for.sens that do not *start* in gr can lead to disjoint
                        # ranges. To be safe, read in training sites covering the union of gr2
                        # and gr (because gr2 can also extend outside of gr). This makes sure all
                        # possible training sites are available AND that only one contigous range
                        # is present in `region`.
                        region=reduce(c(gr, gr2)),
                        integrated.table.path=integrated.table.path,
                        single.cell.id=single.cell.id, quiet=quiet),
                report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            if (nrow(training.sites) == 0) {
                pc <- paste('skipping compute.ab', i, '(empty)')
                ret <- NULL
            } else {
                pc <- perfcheck(paste('compute.ab',i),
                        ab <- compute.ab.given.sites.and.training.data(
                            sites=data.table(chr=as.character(seqnames(gr2)), pos=(end(gr2)+start(gr2)) / 2),
                            training.hsnps=training.sites,
                            ab.fits=ab.fits, quiet=TRUE),  # have to set quiet=TRUE or progress bar will get overridden
                    report.mem=report.mem)
                ret <- cbind(data.table(chr=as.character(seqnames(gr2)), start=start(gr2), end=end(gr2)), ab)
            }
            p(class='sticky', amount=1, pc)

            ret
        })
    }, enable=TRUE)

    return(rbindlist(xs))
}


model.somatic.sensitivity <- function(tiles, muttype=c('snv', 'indel'), alleletype=c('maj', 'min'), sex.chroms=c(), random.seed=0) {
    muttype <- match.arg(muttype)
    alleletype <- match.arg(alleletype)

    # only use the columns corresponding to this muttype/alleletype
    tstr <- paste0(muttype, '.n.training.', alleletype)
    pstr <- paste0(muttype, '.n.training.passed.', alleletype)
    bsc <- paste0('bases.gt.', muttype, '.sc.min.dp')
    bblk <- paste0('bases.gt.', muttype, '.bulk.min.dp')
    nbr <- paste0(muttype, '.n.training.neighborhood')

    model.lhs <- paste0('cbind(', pstr, ', ', tstr, ' - ', pstr, ')')

    # experimental: add log10(depth) to model the negative effect of extreme depth.
    # reasonable depth values is positively related with sensitivity.  now seems to be better
    model.rhs <- paste(
        c('abs(gp.mu)', bsc, bblk, nbr, 'gp.sd', 'mean.sc.dp', 
          #'I(log10(1+mean.sc.dp))',
          'norm.mean.sc.dp',
          # recycle0 - if length(sex.chroms)=0, then paste(.) returns a 0-length result
          paste0('is.', sex.chroms, recycle0=TRUE)),
        collapse=" + ")

    # Build model
    set.seed(random.seed)  # glm may not be deterministic. Just be safe.
    model <- glm(paste(model.lhs, '~', model.rhs), family=binomial, data=tiles[hold.out == FALSE])

    # Predict on the entire data table (INCLUDING training data) and save predictions
    # to caller's table by reference.
    # To ASSESS the model, only consider the predictions on the hold.out==TRUE set.
    predname <- paste0('pred.', muttype, '.', alleletype)
    # 1/(1+exp(-x)) is inverse logit to transform -inf..+inf values to probabilities in [0,1]
    # rankdeficient: important for memory usage. predict() often fails on, e.g., 0-copy
    #       sex chromosomes. when it fails, the default behavior is to attach an attr()
    #       with a named vector of all sites that fail.  for just chrY, this vector is
    #       ~10 Mb in size.  this large attribute is copied to any table that uses the
    #       predname column, so the several copies of the table produced for alleletype=maj,
    #       min and both all receive a copy.  Those copies are further propagated through
    #       glm() calls (which save a copy their input in the model object).
    tiles[, (predname) := 1/(1+exp(-predict(object=model, newdata=.SD, rankdeficient='NA')))]

    model
}


# convenience function to pare down @spatial.sensitivity to columns relevant
# only to a specific (muttype, alleletype) combo. the returned data.table will
# have columns named "pred", "npass", "ntotal" regardless of which (muttype,
# alleletype) pair was passed.
#
# alleletype='both' takes the average predicted sensitivity of major and minor
#   alleles, representing the assumption that somatic mutations are equally
#   likely to occur on either allele. this assumption is not always a good one:
#   e.g., balancer chromosomes, perhaps Barr bodies, etc.
select.muttype.and.alleletype <- function(object, data=object@spatial.sensitivity$data, muttype=c('snv', 'indel'), alleletype=c('both', 'maj', 'min')) {
    muttype <- match.arg(muttype)
    alleletype <- match.arg(alleletype)
    
    if (alleletype == 'both') {
        pred.maj <- paste0('pred.', muttype, '.maj')
        pred.min <- paste0('pred.', muttype, '.min')
        ret <- data[, .(hold.out,
            pred=(get(pred.maj) + get(pred.min)) / 2,
            npass=get(paste0(muttype, '.n.training.passed')),
            ntotal=get(paste0(muttype, '.n.training')),
            ncalls=get(paste0(muttype, '.n.calls')))]
    } else {
        ret <- data[, .(hold.out,
            pred=get(paste0('pred.', muttype, '.', alleletype)),
            npass=get(paste0(muttype, '.n.training.passed.', alleletype)),
            ntotal=get(paste0(muttype, '.n.training.', alleletype)),
            ncalls=get(paste0(muttype, '.n.calls.', alleletype)))]
    }
    ret
}


assess.predicted.somatic.sensitivity <- function(object, muttype=c('snv', 'indel'), alleletype=c('both', 'maj', 'min'))
{
    # For assessing the model, need to use only the held-out tiles that were not
    # used for training.
    # ntotal > 0 makes the comparison similar to genome-wide mean germline sens.
    stratify.tiles.by.somatic.sensitivity(
        data=object@spatial.sensitivity$data[hold.out == TRUE],
        muttype=muttype, alleletype=alleletype)
}


# Return a table for the named covariate with:
#   - binned covariate value
#   - sensitivity for all combinations of:
#       * muttype=snv, indel 
#       * allele=maj, min
#       * sensitivity=measured (at hSNPs), predicted (by model)
#   - number of tiles in genome with this binned cov value
assess.covariate <- function(object, cov, cov.bin.digits=2) {
    xform <- identity
    cov.name <- cov
    if (cov == 'abs(gp.mu)' | cov == 'gp.sd') {
        # abs(gp.mu) must be mapped to gp.mu because the column in the table is named gp.mu
        # the "abs(gp.mu)" name is kept as the name of the covariate (in cov.name) because
        # that's the coef name in the model.
        if (cov == 'abs(gp.mu)')
            cov <- 'gp.mu'

        # abs() does nothing for gp.sd, which is >0. Just want round(2)
        xform <- function(x) round(abs(x),cov.bin.digits)
    } else if (cov == 'mean.sc.dp') {
        xform <- function(x) round(x,0)   # integerize, but nearest
    } else if (cov == 'norm.mean.sc.dp') {
        xform <- function(x) round(x, cov.bin.digits)
    } else if (cov == 'I(log10(1 + mean.sc.dp))') {
        cov <- 'mean.sc.dp'
        xform <- function(x) round(log10(1+x), cov.bin.digits)
    # E.g., when "X" is the name of chromosome X, the model will have a covariate
    # named is.XTRUE.  Code here avoids assuming the name of the sex chromosomes.
    } else if (substr(cov, 1, 3) == 'is.' & substr(cov, nchar(cov)-3, nchar(cov)) == 'TRUE') {
        cov <- 'chr'
        xform <- function(x) x == sub('^is.', '', sub('TRUE$', '', cov.name))
    }

    # is.na(mean.bulk.dp) is essentially an alias for tiles in unassembled genome
    # regions. there are ~71 tiles containing hSNPs with is.na(mean.bulk.dp) vs.
    # 195,046 tiles with no hSNPs in the neighborhood of +/- 10kb.
    object@spatial.sensitivity$data[!is.na(mean.bulk.dp),
        .(cov.name=..cov.name,
          snv.sens.maj=sum(snv.n.training.passed.maj, na.rm=TRUE) / sum(snv.n.training.maj, na.rm=TRUE),
          snv.sens.min=sum(snv.n.training.passed.min, na.rm=TRUE) / sum(snv.n.training.min, na.rm=TRUE),
          pred.snv.sens.maj=mean(pred.snv.maj, na.rm=TRUE),
          pred.snv.sens.min=mean(pred.snv.min, na.rm=TRUE),
          indel.sens.maj=sum(indel.n.training.passed.maj, na.rm=TRUE) / sum(indel.n.training.maj, na.rm=TRUE),
          indel.sens.min=sum(indel.n.training.passed.min, na.rm=TRUE) / sum(indel.n.training.min, na.rm=TRUE),
          pred.indel.sens.maj=mean(pred.indel.maj, na.rm=TRUE),
          pred.indel.sens.min=mean(pred.indel.min, na.rm=TRUE),
          n.tiles=nrow(.SD)),
        by=.(cov=xform(get(cov)))][order(cov)]
}


# predicted sensitivity is rounded to `stratify.sensitivity.digits` decimals. 2 is
# good in practice. increasing the digits will use a closer approximation of the
# model's sensitivity, which could increase accuracy in some cases. however, more
# digits will also divide the genome into more, finer regions, thereby decreasing
# the number of somatic calls and germline sites in each region. this reduction in
# sites increases noise in sensitivity estimates, which may outweigh the increase
# in accuracy gained by using more digits from the predicted sensitivity.
stratify.tiles.by.somatic.sensitivity <- function(object, data=object@spatial.sensitivity$data,
    muttype=c('snv', 'indel'), alleletype=c('both', 'maj', 'min'), stratify.sensitivity.digits=2)
{
    ret <- select.muttype.and.alleletype(data=data, muttype=muttype, alleletype=alleletype)

    ret[, .(npass=sum(npass), ntotal=sum(ntotal),
            sens=sum(npass)/sum(ntotal),
            ncalls=sum(ncalls),
            n=nrow(.SD)), by=.(pred=round(pred, stratify.sensitivity.digits))][order(pred)]
}


# as above, when comparing spatial and germline L-O-O estimates of sensitivity,
# windows with 0 germline sites must be excluded. However, do not need to
# consider hold.out vs. training.
# keep.germline.zeros=TRUE may provide a better genome-wide extrapolation of
# sensitivity because germline L-O-O estimates cannot, by definition, provide
# information about regions of the genome where no germline variants are
# present.
somatic.sensitivity <- function(object, data=object@spatial.sensitivity$data,
    muttype=c('snv', 'indel'), alleletype=c('both', 'maj', 'min'), keep.germline.zeros=FALSE)
{
    if (!missing(object) & !missing(data))
        stop('exactly one of `object` or `data` may be specified')

    muttype <- match.arg(muttype, several.ok=TRUE)
    alleletype <- match.arg(alleletype, several.ok=TRUE)

    setNames(lapply(muttype, function(mt) {
        ret <- do.call(cbind, lapply(alleletype, function(at) {
            this.ret <- select.muttype.and.alleletype(data=data, muttype=mt, alleletype=at)
            unlist(this.ret[!is.na(pred) & (keep.germline.zeros | ntotal > 0),
                .(germline.training=sum(npass)/sum(ntotal),
                  predicted.somatic.mean=mean(pred, na.rm=TRUE),
                  predicted.somatic.stdev=sd(pred, na.rm=TRUE),
                  predicted.somatic.q25=quantile(pred, prob=0.25, na.rm=TRUE),
                  predicted.somatic.q75=quantile(pred, prob=0.75, na.rm=TRUE))
            ])
        }))
        colnames(ret) <- alleletype
        ret
    }), muttype)
}


# *********** EXPERIMENTAL TOTAL MUTATION BURDEN ESTIMATORS ************
# Use an EQUAL-BURDEN ASSUMPTION across major and minor alleles to estimate
# total somatic burden. That is, assume that the same number of mutations
# occur on the major allele (VAF >= 50%) as the minor allele. This assumption
# is helpful because sensitivity on the minor allele is much lower than the
# major allele and leads to noisier estimates. The allele that is major or
# minor changes with genome location; it does not refer to either the maternal
# or paternal allele (in humans).
#
# MAJOR RESERVATIONS: in low mutation burden settings, the regression on
# sensitivity for the negative binomial and poisson models below can fail
# if the few calls happen to land in low sensitivity bins.  Further,
# sensitivity should not have an exponential relationship with the call
# rate as the count models both imply; it should have a simple linear
# relationship.  As such, model (4) might be the best estimate by directly
# estimating the mean mutation rate per tile, assuming all tiles have the
# same mutation rate (which is not true). The main drawback of the mean
# estimator is that it is often skewed by calls in low sensitivity
# bins (which could be false positive calls or poorly predicted sens.) and/or
# small bins.  The other methods include the bin size in the probablistic
# model, which mitigates the extra noise of small bins to some degree.
#
# Finally, even though we try to ensure a model fit succeeds, there remain
# corner cases where one of these models will fail to fit for numerical
# reasons. This is handled by tryCatch() and, when failure occurs, the
# model will be NULL.
#
# 4 methods are used to produce a total burden estimate (i.e., corrected for
# somatic mutation detection sensitivity) by first stratifying the genome into
# (non-contiguous) bins with similar somatic sensitivity ("similar" is controlled
# by stratify.sensitivity.digits). 
#
# In order of worst to best, the methods are:
#   1. Negative binomial simulation. In each bin i, simulate the number of false
#      negative mutations via
#         N_i = # of somatic mutation calls in bin i
#         X_i ~ NegBin(size=N_i, prob=sensitivity_i).
#      The total number of mutations per bin is then X_i + N_i and summing over i
#      gives the total number of mutations in the genome.  This method suffers
#      from the fact that NegBin is not applicable to N_i=0 or sensitivity_i=0,
#      so many bins must be assigned X_i=0.
sim.negbin.estimator <- function(strat.tab, n.sims) {
    # use a negative binomial process to estimate total number of mutations in
    # each sensitivity bin (i.e., number called and number missed(=negbin value))
    # This is usually not a great estimate because of the two conditions below:
    # sens > 0 & ncalls > 0. The areas of the genome with no calls or no sensitivity
    # must be ignored because the negative binomial distribution does not support
    # size=0 or prob=0.
    sims <- rowSums(apply(strat.tab[sens > 0 & ntotal > 0,.(ncalls, sens)], 1, function(row)
        if (row[1] == 0) {
            rep(0, n.sims)
        } else {
            row[1] + rnbinom(n=n.sims, size=row[1], prob=row[2])
        }
    ))

    list(
        n.sims=n.sims,
        burden=mean(sims),
        mean=mean(sims), median=median(sims),
        ci95=quantile(sims, probs=c(0.025, 0.975)),
        ci.99=quantile(sims, probs=c(0.005, 0.995)),
        sims=sims)
}
#   2. Poisson regression on N_i / S_i ~ sensitivity_i, where S_i is the number
#      of tiles in bin i.  That is, the rate of mutations per tile is being
#      modeled.  The primary assumption is that spatial sensitivity is independent
#      of any spatial bias in mutation distribution.  E.g., in cancer mutations
#      are known to be spatial biased toward (i.e., overrepresented) in closed
#      chromatin regions.
poisson.estimator <- function(strat.tab) {
    pois.m <- stats::glm(ncalls ~ sens + offset(log(n)),
        # ntotal > 0: for low depth cells, can have cases where 0 germline sites
        # are within bins with a similar sensitivity prediction.  Leads to sens=NaN.
        data=strat.tab[sens > 0 & ntotal > 0],
        family=poisson)

    list(model=pois.m,
        burden=exp(predict(object=pois.m, newdata=data.frame(sens=1, n=sum(strat.tab$n)))))
        # The total extrapolation uses a bin the size of the whole genome (sum(strat.tab$n))
        # with detection sensitivity = 100%.  The rate parameter estimated by the model is
}
#   3. Negative binomial regression on N_i / S_i ~ sensitivity_i.  This is the
#      same model as (2), but with an added overdispersion parameter.
negbin.estimator <- function(strat.tab) {
    data <- strat.tab[sens > 0 & ntotal > 0]
    # glm.nb fails if all of ncalls=0. glm(family=poisson) does not have this issue
    if (sum(data$ncalls) == 0) {
        negbin.m <- NULL
        burden <- 0
    } else {
        negbin.m <- tryCatch(MASS::glm.nb(ncalls ~ sens + offset(log(n)), data=data),
            error=function(e) { cat(paste("ERROR occurred in MASS::glm.nb (see below), skipping and returning NULL model\n", e) ); return(NULL) })
        burden <- 0
        if (!is.null(negbin.m))
            burden <- exp(predict(object=negbin.m, newdata=data.frame(sens=1, n=sum(strat.tab$n))))
    }
    list(model=negbin.m, burden=burden)
}
#   4. Mean estimator.
mean.estimator <- function(strat.tab) {
    # Don't use pred=NA, pred=0 or sens=0 bins to estimate the per-tile rate,
    # but do use all bins in the extrapolation (sum(strat.tab$n)).
    rate.per.tile.per.sens.bin <- strat.tab[!is.na(pred) & pred > 0 & sens > 0 & ntotal > 0, ncalls / sens / n]
    rate.per.tile <- mean(rate.per.tile.per.sens.bin)
    list(
        rate.per.tile.per.sens.bin=rate.per.tile.per.sens.bin,
        rate.per.tile=rate.per.tile,
        burden=rate.per.tile*sum(strat.tab$n)
    )
}
#
# N.B. in this context 'major' and 'minor' alleles refer amplification bias:
# the allele that was more amplified as "major" and less amplified is "minor".
# there is no differentiation between VAF>=50% due to sampling noise and
# actual differences in amplification level.
estimate.burden.by.spatial.sensitivity <- function(object, data=object@spatial.sensitivity$data,
     muttype=c('snv', 'indel'), n.sims=1e4, stratify.sensitivity.digits=2, random.seed=0)
{
    set.seed(random.seed)
    muttype <- match.arg(muttype)

    # use the model-assigned sensitivity to every tile in the genome to roughly
    # group the genome by somatic calling sensitivity. once the model's predicted
    # sensitivity is used to create an aggregated region for each sensitivity value,
    # the GERMLINE sites in each region are used to calculate sensitivity.
    allele.types <- c('both', 'maj', 'min')
    strat.tabs <- setNames(lapply(allele.types, function(at) {
        strat.tab <- stratify.tiles.by.somatic.sensitivity(
            data=data, muttype=muttype, alleletype=at,
            stratify.sensitivity.digits=stratify.sensitivity.digits)
    }), allele.types)

    # Estimate 1.
    sim.negbin.estimates <- lapply(strat.tabs, sim.negbin.estimator, n.sims=n.sims)
    # Estimate 2
    poisson.estimates <- lapply(strat.tabs, poisson.estimator)
    # Estimate 3
    negbin.estimates <- lapply(strat.tabs, negbin.estimator)
    # Estimate 4
    mean.estimates <- lapply(strat.tabs, mean.estimator)

    # EQUAL-BURDEN ASSUMPTION: with a decent estimate of the number of mutations
    # on the VAF>=50% (major) allele, assume the other allele (diploid assumption,
    # obviously) has the same number of mutations.
    sim.negbin.estimates$equal.burden <- sim.negbin.estimates$maj$burden*2
    poisson.estimates$equal.burden <- poisson.estimates$maj$burden*2
    negbin.estimates$equal.burden <- negbin.estimates$maj$burden*2
    mean.estimates$equal.burden <- mean.estimates$maj$burden*2

    # return the simulations, metadata and some useful precomputed summaries
    list(muttype=muttype, random.seed=random.seed,
        strat.tabs=strat.tabs,
        sim.negbin.estimates=sim.negbin.estimates,
        poisson.estimates=poisson.estimates,
        negbin.estimates=negbin.estimates,
        mean.estimates=mean.estimates)
}


experimental.mutburden.estimators <- function(object, muttype, estimators=c('sim.negbin.estimates', 'poisson.estimates', 'negbin.estimates', 'mean.estimates')) {
    sapply(object@spatial.sensitivity$burden[[muttype]][estimators], function(x)
        c(maj=x$maj$burden, min=x$min$burden, both=x$both$burden, equal.assumption=x$equal.burden))
}


compare.mutburden.estimators <- function(object, muttype=c('snv', 'indel')) {
    muttype <- match.arg(muttype, several.ok=TRUE)

    setNames(lapply(muttype, function(mt) {
        original.estimator <- c(maj=NA, min=NA, both=mutburden(object, mt)[[1]], equal.assumption=NA)
        spatial.estimators <- experimental.mutburden.estimators(object, mt)
        round(rbind(original.estimator, t(spatial.estimators)), 1)
    }), muttype)
}


# sites should be filtered for training.site==TRUE and one muttype
# N.B. allele balance is probably a better indicator of underlying
# amplification imbalance than VAF.
count.germline.sites.for.sens <- function(grs, sites, seqinfo, neighborhood.tiles) {
    sites.to.gr <- function(sites) {
        GenomicRanges::GRanges(seqnames=sites$chr,
                ranges=IRanges::IRanges(start=sites$pos, width=1),
                seqinfo=seqinfo)
    }
    ret <- data.table::data.table(
        n.training=IRanges::countOverlaps(grs,
            sites.to.gr(sites)),
        n.training.maj=IRanges::countOverlaps(grs,
            sites.to.gr(sites[af >= 0.5])),
        n.training.min=IRanges::countOverlaps(grs,
            sites.to.gr(sites[af < 0.5])),
        n.training.passed=IRanges::countOverlaps(grs,
            sites.to.gr(sites[training.pass == TRUE])),
        n.training.passed.maj=IRanges::countOverlaps(grs,
            sites.to.gr(sites[training.pass == TRUE & af >= 0.5])),
        n.training.passed.min=IRanges::countOverlaps(grs,
            sites.to.gr(sites[training.pass == TRUE & af < 0.5])))
    ret[, n.training.neighborhood :=
        as.integer(stats::filter(n.training, filter=rep(1, 2*neighborhood.tiles), sides=2))]
    ret
}


# should really reuse above code. oh well.
count.somatic.sites.for.sens <- function(grs, sites, seqinfo, neighborhood.tiles) {
    sites.to.gr <- function(sites) {
        GenomicRanges::GRanges(seqnames=sites$chr,
                ranges=IRanges::IRanges(start=sites$pos, width=1),
                seqinfo=seqinfo)
    }
    ret <- data.table::data.table(
        n.calls=IRanges::countOverlaps(grs, sites.to.gr(sites)),
        n.calls.maj=IRanges::countOverlaps(grs, sites.to.gr(sites[af >= 0.5])),
        n.calls.min=IRanges::countOverlaps(grs, sites.to.gr(sites[af < 0.5])))
    ret[, n.neighborhood :=
        as.integer(stats::filter(n.calls, filter=rep(1, 2*neighborhood.tiles), sides=2))]
    ret
}

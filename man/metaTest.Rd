\name{metaTest}
\alias{metaTest}
\title{Meta-analysis using several RNA-Seq statistics}
\usage{
    metaTest(cpList,
        metaP = c("simes", "bonferroni", "fisher", "harmonic",
        "dperm_min", "dperm_max", "dperm_weight", "fperm", 
        "whitlock", "minp", "maxp", "weight", "pandora", 
        "none"), counts, sampleList, statistics, statArgs, 
        libsizeList, nperm = 10000, 
        weight = rep(1/length(statistics), length(statistics)), 
        pOffset = NULL, rc = NULL)
}
\arguments{
    \item{cpList}{a named list whose names are the contrasts
    requested from metaseqr2. Each member is a p-value matrix
    whose colnames are the names of the statistical tests
    applied to the data. See the main \code{\link{metaseqr2}}
    help page.}

    \item{metaP}{the p-value combination method to use. See
    the main \code{\link{metaseqr2}} help page.}

    \item{counts}{the normalized and possibly filtered read
    counts matrix. See the main \code{\link{metaseqr2}} help
    page.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition. See the main
    \code{\link{metaseqr2}} help page.}

    \item{statistics}{the statistical algorithms used in
    metaseqr2. See the main \code{\link{metaseqr2}} help page.}

    \item{statArgs}{the parameters for each statistical
    argument. See the main \code{\link{metaseqr2}} help page.}

    \item{libsizeList}{a list with library sizes. See the
    main \code{\link{metaseqr2}} and the \code{stat*} help
    pages.}

    \item{nperm}{the number of permutations (Monte Carlo
    simulations) to perform.}

    \item{weight}{a numeric vector of weights for each
    statistical algorithm.}
    
    \item{pOffset}{\code{NULL} (default) or a fixed 
    numeric value between 0 and 1. See also the main 
    \code{\link{metaseqr2}} man page.}

    \item{rc}{the fraction of the available cores to use
    in a multicore system.}
}
\value{
    A named list with combined p-values. The names are the
    contrasts and the list members are combined p-value
    vectors, one for each contrast.
}
\description{
    This function calculates the combined p-values when
    multiple statistical algorithms are applied to the input
    dataset. It is a helper and it requires very specific
    arguments so it should not be used individually
}
\details{
    Ideally one would want to create the same set of indices 
    for a given dataset so as to create reproducible p-values. 
    To achieve this, use the \code{set.seed} function prior 
    to any calculations.
}
\examples{
cpList <- list(a=matrix(runif(100),50,2))
metaP <- metaTest(cpList,"simes")
}
\author{
    Panagiotis Moulos
}

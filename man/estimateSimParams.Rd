\name{estimateSimParams}
\alias{estimateSimParams}
\title{Estimate negative binomial parameters from 
    real data}
\usage{
    estimateSimParams(realCounts, libsizeGt = 3e+6,
        rowmeansGt = 5,eps = 1e-11, rc = NULL, draw = FALSE)
}
\arguments{
    \item{realCounts}{a text tab-delimited file 
    with real RNA-Seq data. See Details.}

    \item{libsizeGt}{a library size below which 
    samples are excluded from parameter estimation 
    (default: \code{3000000}).}

    \item{rowmeansGt}{a row means (mean counts 
    over samples for each gene) below which 
    genes are excluded from parameter estimation
    (default: 5).}

    \item{eps}{the tolerance for the convergence 
    of \code{\link{optimize}} function. Defaults 
    to 1e-11.}

    \item{rc}{in case of parallel optimization, the 
    fraction of the available cores to use.}

    \item{draw}{boolean to determine whether to 
    plot the estimated simulation parameters 
    (mean and dispersion) or not. Defaults to 
    \code{FALSE} (do not draw a mean-dispersion 
    scatterplot).}
}
\value{
    A named list with two members: \code{muHat}
    which contains negative binomial mean 
    estimates and \code{phiHat} which contains 
    dispersion estimates.
}
\description{
    This function reads a read counts table 
    containing real RNA-Seq data (preferebly 
    with more than 20 samples so as to get as 
    much accurate as possible estimations) and 
    calculates a population of count means and 
    dispersion parameters which can be used to 
    simulate an RNA-Seq dataset with synthetic 
    genes by drawing from a negative binomial 
    distribution. This function works in the 
    same way as described in (Soneson and 
    Delorenzi, BMC Bioinformatics, 2013) and 
    (Robles et al., BMC Genomics, 2012).
}
\details{
    Regarding \code{realCounts}, the file should strictly
    contain a unique gene name (e.g. Ensembl accession) in 
    the first column and all other columns should contain 
    read counts for each gene. Each column must be named
    with a unique sample identifier. See examples in the 
    ReCount database 
    \url{http://bowtie-bio.sourceforge.net/recount/}.
    
    Also, the parameter estimation involves a lot of 
    random sampling. For guaranteed reproducibility,
    be sure to use \code{set.seed} prior to any
    calculations. By default, when the metaseqR2 package
    is loaded, the seed is set to \code{42}.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
parList <- estimateSimParams(dataMatrix,libsizeGt=3e+4)
}
\author{
    Panagiotis Moulos
}


\name{makeSimDataSd}
\alias{makeSimDataSd}
\title{Create simulated counts using the 
    Soneson-Delorenzi method}
\usage{
    makeSimDataSd(N, param, samples = c(5, 5),
        ndeg = rep(round(0.1*N), 2), fcBasis = 1.5,
        libsizeRange = c(0.7, 1.4), libsizeMag = 1e+7,
        modelOrg = NULL, simLengthBias = FALSE, 
        seed = NULL)
}
\arguments{
    \item{N}{the number of genes to produce.}

    \item{param}{a named list with negative binomial 
    parameter sets to sample from. The first member is
    the mean parameter to sample from (\code{muHat}) 
    and the second the dispersion (\code{phiHat}). 
    This list can be created with the 
    \code{\link{estimateSimParams}} function.}

    \item{samples}{a vector with 2 integers, 
    which are the number of samples for each 
    condition (two conditions currently supported).}

    \item{ndeg}{a vector with 2 integers, which are 
    the number of differentially expressed genes to 
    be produced. The first element is the number of 
    up-regulated genes while the second is the 
    number of down-regulated genes.}

    \item{fcBasis}{the minimum fold-change for 
    deregulation.}

    \item{libsizeRange}{a vector with 2 numbers 
    (generally small, see the default), as they 
    are multiplied with \code{libsizeMag}. These 
    numbers control the library sized of the 
    synthetic data to be produced.}

    \item{libsizeMag}{a (big) number to multiply 
    the \code{libsizeRange} to produce library 
    sizes.}

    \item{modelOrg}{the organism from which the 
    real data are derived from. It must be one 
    of the supported organisms (see the main 
    \code{\link{metaseqr2}} help page). It is used 
    to sample real values for GC content.}

    \item{simLengthBias}{a boolean to instruct 
    the simulator to create genes whose read counts is
    proportional to their length. This is achieved by 
    sorting in increasing order the mean parameter of 
    the negative binomial distribution (and the 
    dispersion according to the mean) which will cause 
    an increasing gene count length with the sampling. 
    The sampled lengths are also sorted so that in the 
    final gene list, shorter genes have less counts as 
    compared to the longer ones. The default is FALSE.}

    \item{seed}{a seed to use with random number 
    generation for reproducibility.}
}
\value{
    A named list with two members. The first 
    member (\code{simdata}) contains the 
    synthetic dataset 
}
\description{
    This function creates simulated RNA-Seq gene 
    expression datasets using the method presented 
    in (Soneson and Delorenzi, BMC Bioinformatics, 
    2013). For the time being, it creates only 
    simulated datasets with two conditions.
}
\examples{
\donttest{
require(DESeq)
dataMatrix <- counts(makeExampleCountDataSet())
## File "bottomly_read_counts.txt" from the ReCount database
#download.file(paste("http://bowtie-bio.sourceforge.net/recount/",
#    "countTables/bottomly_count_table.txt",sep=""),
#    destfile="~/bottomly_count_table.txt")
N <- 10000
#parList <- estimateSimParams("~/bottomly_read_counts.txt")
parList <- estimateSimParams(dataMatrix)
sim <- makeSimDataSd(N,parList)
synthData <- sim$simdata
trueDeg <- which(sim$truedeg!=0)
}
}
\author{
    Panagiotis Moulos
}


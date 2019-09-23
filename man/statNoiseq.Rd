\name{statNoiseq}
\alias{statNoiseq}
\title{Statistical testing with NOISeq}
\usage{
    statNoiseq(object, sampleList, contrastList = NULL,
        statArgs = NULL, geneData = NULL, logOffset = 1)
}
\arguments{
    \item{object}{a matrix or an object specific to each
    normalization algorithm supported by metaseqR2, containing
    normalized counts. See also Details.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{contrastList}{vector of contrasts as defined in the 
    main help page of \code{\link{metaseqr2}}. See also 
    Details.}

    \item{statArgs}{a list of edgeR statistical algorithm
    parameters. See the result of
    \code{getDefaults("statistics",} \code{"noiseq")} for an
    example and how you can modify it.}

    \item{geneData}{an optional annotation data frame (such
    the ones produced by \code{get.annotation} which contains
    the GC content for each gene and from which the gene
    lengths can be inferred by chromosome coordinates.}

    \item{logOffset}{a number to be added to each element of
    data matrix in order to avoid Infinity on log type data
    transformations.}
}
\value{
    A named list of NOISeq q-values, whose names are the
    names of the contrasts.
}
\description{
    This function is a wrapper over NOISeq statistical
    testing. It accepts a matrix of normalized gene counts or
    an S4 object specific to each normalization algorithm
    supported by metaseqR2.
}
\details{
    Regarding \code{object}, apart from \code{matrix} (also 
    for NOISeq), the object can be a \code{SeqExpressionSet} 
    (EDASeq), \code{CountDataSet} (DESeq), \code{DGEList} 
    (edgeR), \code{DESeqDataSet} (DESeq2), \code{SeqCountSet} 
    (DSS) or \code{ABSDataSet} (ABSSeq).
    
    Regarding \code{contrastList} it can also be a named 
    structured list of contrasts as returned by the internal
    function \code{metaseqR2:::makeContrastList}.
}
\examples{
\donttest{
require(DESeq)
dataMatrix <- counts(makeExampleCountDataSet())
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
lengths <- round(1000*runif(nrow(dataMatrix)))
starts <- round(1000*runif(nrow(dataMatrix)))
ends <- starts + lengths
gc=runif(nrow(dataMatrix))
geneData <- data.frame(
    chromosome=c(rep("chr1",nrow(dataMatrix)/2),
        rep("chr2",nrow(dataMatrix)/2)),
    start=starts,end=ends,gene_id=rownames(dataMatrix),gc_content=gc
)
norm.dataMatrix <- normalize.noiseq(dataMatrix,sampleList,geneData)
p <- statNoiseq(norm.dataMatrix,sampleList,contrast,
    geneData=geneData)
}
}
\author{
    Panagiotis Moulos
}


\name{statBayseq}
\alias{statBayseq}
\title{Statistical testing with baySeq}
\usage{
    statBayseq(object, sampleList, contrastList = NULL,
        statArgs = NULL, libsizeList = NULL)
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

    \item{statArgs}{a list of bayseq statistical algorithm
    parameters. See the result of
    \code{getDefaults("statistics",} \code{"bayseq")} for an
    example and how you can modify it.}

    \item{libsizeList}{an optional named list where names
    represent samples (MUST be the same as the samples in
    \code{sampleList}) and members are the library sizes
    (the sequencing depth) for each sample. If not provided,
    they will be estimated from baySeq.}
}
\value{
    A named list of the value 1-likelihood that a gene is
    differentially expressed, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over baySeq statistical
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
# Careful on testing, baySeq is slow
require(DESeq)
require(baySeq)
dataMatrix <- DESeq::counts(makeExampleCountDataSet())
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
normDataMatrix <- normalizeEdger(dataMatrix,sampleList)
p <- statBayseq(normDataMatrix,sampleList,contrast)
}
}
\author{
    Panagiotis Moulos
}


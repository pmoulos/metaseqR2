\name{statEdger}
\alias{statEdger}
\title{Statistical testing with edgeR}
\usage{
    statEdger(object, sampleList, contrastList = NULL,
        statArgs = NULL)
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
    \code{getDefaults("statistics",} \code{"edger")} for an
    example and how you can modify it.}
}
\value{
    A named list of p-values, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over edgeR statistical
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
normDataMatrix <- normalize.edger(dataMatrix,sampleList)
p <- statEdger(normDataMatrix,sampleList,contrast)
}
}
\author{
    Panagiotis Moulos
}


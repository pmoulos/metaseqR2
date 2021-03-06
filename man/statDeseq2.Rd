\name{statDeseq2}
\alias{statDeseq2}
\title{Statistical testing with DESeq2}
\usage{
    statDeseq2(object, sampleList, contrastList = NULL,
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

    \item{statArgs}{a list of DESeq statistical algorithm
    parameters. See the result of
    \code{getDefaults("statistics",} \code{"deseq")} for an
    example and how you can modify it. It is not required
    when the input object is already a CountDataSet from
    DESeq normalization as the dispersions are already
    estimated.}
}
\value{
    A named list of p-values, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over DESeq statistical
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
require(DESeq2)
dataMatrix <- metaseqR2:::exampleCountData(1000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
normDataMatrix <- normalizeDeseq2(dataMatrix,sampleList)
p <- statDeseq2(normDataMatrix,sampleList,contrast)
}
\author{
    Dionysios Fanidis
}


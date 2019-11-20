\name{statNbpseq}
\alias{statNbpseq}
\title{Statistical testing with NBPSeq}
\usage{
    statNbpseq(object, sampleList, contrastList = NULL,
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

    \item{statArgs}{a list of NBPSeq statistical algorithm
    parameters. See the result of
    \code{getDefaults("statistics",} \code{"nbpseq")} for an
    example and how you can modify it. It is not required
    when the input object is already a list from NBPSeq
    normalization as the dispersions are already estimated.}

    \item{libsizeList}{an optional named list where names
    represent samples (MUST be the same as the samples
    \code{in sampleList}) and members are the library sizes
    (the sequencing depth) for each sample. If not provided,
    the default is the column sums of the \code{object}
    matrix.}
}
\value{
    A named list of p-values, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over NBPSeq statistical
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
\note{
    There is currently a problem with the NBPSeq package and
    the workflow that is specific to the NBPSeq package. The
    problem has to do with function exporting as there are
    certain functions which are not recognized from the
    package internally. For this reason and until it is
    fixed, only the Smyth workflow will be available with the
    NBPSeq package.
}
\examples{
require(DESeq)
require(NBPSeq)
dataMatrix <- DESeq::counts(makeExampleCountDataSet())
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
normDataMatrix <- normalizeNbpseq(dataMatrix,sampleList)
p <- statNbpseq(normDataMatrix,sampleList,contrast)
}
\author{
    Panagiotis Moulos
}


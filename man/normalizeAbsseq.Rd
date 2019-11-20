\name{normalizeAbsseq}
\alias{normalizeAbsseq}
\title{Normalization based on the ABSSeq package}
\usage{
    normalizeAbsseq(geneCounts, sampleList,
        normArgs = NULL, output = c("matrix", "native"))
}
\arguments{
    \item{geneCounts}{a table where each row represents a
    gene and each column a sample. Each cell contains the
    read counts for each gene and sample. Such a table can be
    produced outside metaseqr2 and is imported during the
    basic metaseqr2 workflow.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{normArgs}{a list of DESeq normalization
    parameters. See the result of
    \code{getDefaults("normalization",} \code{"deseq")} for
    an example and how you can modify it.}

    \item{output}{the class of the output object. It can be
    \code{"matrix"} (default) for versatility with other
    tools or \code{"native"} for the ABSSeq native S4 object
    (ABSDataSet). In the latter case it should be handled
    with suitable ABSSeq methods.}
}
\value{
    A matrix or a ABSDataSet with normalized counts.
}
\description{
    This function is a wrapper over ABSSeq normalization. It
    accepts a matrix of gene counts (e.g. produced by
    importing an externally generated table of counts to the
    main metaseqr2 pipeline).
}
\examples{
require(DESeq)
dataMatrix <- counts(makeExampleCountDataSet())
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotBoxplot(dataMatrix,sampleList)

normDataMatrix <- normalizeAbsseq(dataMatrix,sampleList)
diagplotBoxplot(normDataMatrix,sampleList)
}
\author{
    Dionysios Fanidis
}

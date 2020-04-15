\name{normalizeDss}
\alias{normalizeDss}
\title{Normalization based on the DSS package}
\usage{
    normalizeDss(geneCounts, sampleList,
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
    tools or \code{"native"} for the DSS native S4 object
    (SeqCountSet). In the latter case it should be handled
    with suitable ABSSeq methods.}
}
\value{
    A matrix or a SeqCountSet with normalized counts.
}
\description{
    This function is a wrapper over ABSSeq normalization. It
    accepts a matrix of gene counts (e.g. produced by
    importing an externally generated table of counts to the
    main metaseqr2 pipeline).
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotBoxplot(dataMatrix,sampleList)

normDataMatrix <- normalizeDss(dataMatrix,sampleList)
diagplotBoxplot(normDataMatrix,sampleList)
}
\author{
    Dionysios Fanidis
}

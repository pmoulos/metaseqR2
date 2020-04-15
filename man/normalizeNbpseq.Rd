\name{normalizeNbpseq}
\alias{normalizeNbpseq}
\title{Normalization based on the NBPSeq package}
\usage{
    normalizeNbpseq(geneCounts, sampleList,
        normArgs = NULL, libsizeList = NULL,
        output = c("matrix", "native"))
}
\arguments{
    \item{geneCounts}{a table where each row represents a
    gene and each column a sample. Each cell contains the
    read counts for each gene and sample. Such a table can be
    produced outside metaseqr2 and is imported during the
    basic metaseqr2 workflow.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{normArgs}{a list of NBPSeq normalization
    parameters. See the result of
    \code{getDefaults("normalization",} \code{"nbpseq")} for
    an example and how you can modify it.}

    \item{libsizeList}{an optional named list where names
    represent samples (MUST be the same as the samples in
    \code{sampleList}) and members are the library sizes
    (the sequencing depth) for each sample. If not provided,
    the default is the column sums of the \code{geneCounts}
    matrix.}

    \item{output}{the class of the output object. It can be
    \code{"matrix"} (default) for versatility with other
    tools or \code{"native"} for the NBPSeq native S4 object
    (a specific list). In the latter case it should be
    handled with suitable NBPSeq methods.}
}
\value{
    A matrix with normalized counts or a list with the
    normalized counts and other NBPSeq specific parameters.
}
\description{
    This function is a wrapper over DESeq normalization. It
    accepts a matrix of gene counts (e.g. produced by
    importing an externally generated table of counts to the
    main metaseqr2 pipeline).
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotBoxplot(dataMatrix,sampleList)

normDataMatrix <- normalizeNbpseq(dataMatrix,sampleList)
diagplotBoxplot(normDataMatrix,sampleList)
}
\author{
    Panagiotis Moulos
}

\name{normalizeEdaseq}
\alias{normalizeEdaseq}
\title{Normalization based on the EDASeq package}
\usage{
    normalizeEdaseq(geneCounts, sampleList,
        normArgs = NULL, geneData = NULL,
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

    \item{normArgs}{a list of EDASeq normalization
    parameters. See the result of
    \code{getDefaults("normalization",} \code{"edaseq")} for
    an example and how you can modify it.}

    \item{geneData}{an optional annotation data frame (such
    the ones produced by \code{getAnnotation}) which
    contains the GC content for each gene and from which the
    gene lengths can be inferred by chromosome coordinates.}

    \item{output}{the class of the output object. It can be
    \code{"matrix"} (default) for versatility with other
    tools or \code{"native"} for the EDASeq native S4 object
    (SeqExpressionSet). In the latter case it should be
    handled with suitable EDASeq methods.}
}
\value{
    A matrix or a SeqExpressionSet with normalized counts.
}
\description{
    This function is a wrapper over EDASeq normalization. It
    accepts a matrix of gene counts (e.g. produced by
    importing an externally generated table of counts to the
    main metaseqr2 pipeline).
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotBoxplot(dataMatrix,sampleList)

lengths <- round(1000*runif(nrow(dataMatrix)))
starts <- round(1000*runif(nrow(dataMatrix)))
ends <- starts + lengths
gc=runif(nrow(dataMatrix))
geneData <- data.frame(
    chromosome=c(rep("chr1",nrow(dataMatrix)/2),
        rep("chr2",nrow(dataMatrix)/2)),
    start=starts,end=ends,gene_id=rownames(dataMatrix),gc_content=gc,
    row.names=rownames(dataMatrix)
)
normDataMatrix <- normalizeEdaseq(dataMatrix,sampleList,
    geneData=geneData)
diagplotBoxplot(normDataMatrix,sampleList)
}
\author{
    Panagiotis Moulos
}

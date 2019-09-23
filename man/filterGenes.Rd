\name{filterGenes}
\alias{filterGenes}
\title{Filter gene expression based on gene counts}
\usage{
    filterGenes(geneCounts, geneData, geneFilters,
    sampleList)
}
\arguments{
    \item{geneCounts}{a matrix of gene counts, preferably
    after the normalization procedure.}

    \item{geneData}{an annotation \code{GenomicRanges} usually
    obtained with \code{\link{getAnnotation}} or
    \code{\link{loadAnnotation}}.}

    \item{geneFilters}{a named list with gene filters and
    their parameters. See the main help page of
    \code{\link{metaseqr2}} for details.}
    
    \item{sampleList}{the list containing condition names 
    and the samples under each condition.}
}
\value{
    a named list with three members. The first member 
    (\code{result} is a named list whose names are the 
    gene filter names and its members are the filtered 
    rownames of \code{geneData}. The second member 
    (\code{cutoff} is a named list whose names are the 
    gene filter names and its members are the cutoff 
    values corresponding to each filter. The third member 
    is a matrix of binary flags (0 for non-filtered, 1 
    for filtered) for each gene. The rownames of the 
    flag matrix correspond to gene ids.
}
\description{
    This function performs the gene expression filtering
    based on gene read counts and a set of gene filter rules.
    For more details see the main help pages of
    \code{\link{metaseqr2}}.
}
\examples{
\donttest{
# FIXME: Make working example
data("mm9.geneData",package="metaseqR")
geneCounts <- mm9.geneCounts
sampleList <- sampleList.mm9
geneCounts <- normalize.edger(as.matrix(geneCounts[,9:12]),
    sampleList)
geneData <- get.annotation("mm9","gene")
geneFilters <- get.defaults("gene.filter","mm9")
filter.results <- filterGenes(geneCounts,geneData,
    geneFilters,sampleList)
}
}
\author{
    Panagiotis Moulos
}

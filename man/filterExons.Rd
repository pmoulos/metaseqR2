\name{filterExons}
\alias{filterExons}
\title{Filter gene expression based on exon counts}
\usage{
    filterExons(theCounts, geneData, sampleList,
        exonFilters, rc = 0.8)
}
\arguments{
    \item{theCounts}{a named list created with the
    \code{\link{constructGeneModel}} function.}

    \item{geneData}{an annotation \code{GenomicRanges} usually
    obtained with \code{\link{getAnnotation}} or
    \code{\link{loadAnnotation}}.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{exonFilters}{a named list with exon filters and
    their parameters. See the main help page of
    \code{\link{metaseqr2}} for details.}

    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).}
}
\value{
    a named list with two members. The first member 
    (\code{result} is a named list whose names are the 
    exon filter names and its members are the filtered 
    rownames of \code{geneData}. The second member is a 
    matrix of binary flags (0 for non-filtered, 1 for 
    filtered) for each gene. The rownames of the flag 
    matrix correspond to gene ids.
}
\description{
    This function performs the gene expression filtering
    based on exon read counts and a set of exon filter rules.
    For more details see the main help pages of
    \code{\link{metaseqr2}}.
}
\examples{
\donttest{
# FIXME: Make working example.
data("hg19.exon.data",package="metaseqR")
exon.counts <- hg19.exon.counts
geneData <- get.annotation("hg19","gene")
sampleList <- sampleList.hg19
exonFilters <- get.defaults("exon.filter")
theCounts <- construct.gene.model(exon.counts,sampleList,
    geneData)
filter.results <- filterExons(theCounts,geneData,
    sampleList,exonFilters)
}
}
\author{
    Panagiotis Moulos
}


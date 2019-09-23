\name{constructGeneModel}
\alias{constructGeneModel}
\title{Assemble a gene model based on exon or 3' UTR counts}
\usage{
    constructGeneModel(countData, annoData, type,
        rc = NULL)
}
\arguments{
    \item{countData}{a read counts data frame produced by
    reading the an external file containing the counts or
    by \code{\link{read2count}}.}

    \item{annoData}{an annotation \code{GenomicRanges} object
    from the same organism as \code{countData} (such the ones 
    produced by \code{\link{loadAnnotation}} and 
    \code{\lin{loadAnnotation}}).}
    
    \item{type}{\code{"exon"} to construct a gene count model
    from exon counts or \code{"utr"} to construct gene count
    model from 3' UTR counts.}

    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).}
}
\value{
    A named list where names represent samples. Each list
    member is also a named list where names correspond to
    gene ids and members are named vectors. Each vector is
    named according to the exons corresponding to each gene
    and contains the read counts for each exon. An attribute
    is also assigned to this list with the sum of lengths of
    each exon/3' UTR (to be used later in normalization 
    factors, e.g. when outputing \code{"rpgm"} values). This
    structure is also used for exon filtering.
}
\description{
    This function assembles gene models (single genes, not
    isoforms) based on the input exon/3' UTR read counts and a
    related annotation \code{GenomicRanges} object. The
    \code{annoData} argument should have a specific format
    and for this reason it's better to use 
    \code{\link{getAnnotation}} or \code{\link{loadAnnotation}}. 
    This function is intended mostly for internal use but can be 
    used if the requirements are met.
}
\examples{
\donttest{
# FIXME: Include an exon counts table and an annotation
#        ranges in examples.
# Takes some time to run...
data("hg19.exon.data",package="metaseqR")
gene.data <- get.annotation("hg19","gene","ensembl")
reduced.gene.data <- reduce.gene.data(hg19.exon.counts,
    gene.data)
multic <- check.parallel(0.4)
gene.model <- constructGeneModel(hg19.exon.counts,
    sample.list.hg19,gene.data,multic)
}
}
\author{
    Panagiotis Moulos
}

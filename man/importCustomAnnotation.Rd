\name{importCustomAnnotation}
\alias{importCustomAnnotation}
\title{Import a metaseqR2 custom annotation element}
\usage{
    importCustomAnnotation(gtfFile, metadata,
        level = c("gene", "transcript", "exon"),
	    type = c("gene", "exon", "utr"))
}
\arguments{
    \item{gtfFile}{a GTF file containing the gene structure
    of the organism to be imported.}

    \item{metadata}{a list with additional information about
    the annotation to be imported. The same as in the
    \code{\link{buildCustomAnnotation}} man page.}
    
    \item{level}{same as the \code{transLevel} in 
    \code{\link{metaseqr2}}.}
    
    \item{type}{same as the \code{countType} in 
    \code{\link{metaseqr2}}.}
    
    \item{version}{same as the \code{version} in 
    \code{\link{metaseqr2}}.}
    
    \item{db}{same as the \code{db} in 
    \code{\link{buildAnnotationDatabase}}.}
    
    \item{rc}{same as the \code{rc} in 
    \code{\link{buildAnnotationDatabase}}.}
}
\value{
    The function returns a \code{GenomicRanges} object with
    the requested annotation.
}
\description{
    This function creates a local annotation database to be
    used with metaseqr2 so as to avoid long time on the fly 
    annotation downloads and formatting.
}
\examples{
    # FIXME: Example
}
\author{
    Panagiotis Moulos
}
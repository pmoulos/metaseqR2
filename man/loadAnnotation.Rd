\name{loadAnnotation}
\alias{loadAnnotation}
\title{Load a metaseqR2 annotation element}
\usage{
    loadAnnotation(genome, refdb, 
        level = c("gene", "transcript", "exon"),
	    type = c("gene", "exon", "utr"), version="auto",
	    db = file.path(system.file(package = "metaseqR2"),
	        "annotation.sqlite"),summarized = FALSE, 
	        asdf = FALSE, rc = NULL)
}
\arguments{
    \item{genome}{a \code{\link{metaseqr2}} supported
    organisms or a custom, imported by the user, name. See 
    also the main \code{\link{metaseqr2}} man page.}

    \item{refdb}{a \code{\link{metaseqr2}} supported
    annotation source or a custom, imported by the user, name.
    See also the main \code{\link{metaseqr2}} man page.}
    
    \item{level}{same as the \code{transLevel} in 
    \code{\link{metaseqr2}}.}
    
    \item{type}{same as the \code{countType} in 
    \code{\link{metaseqr2}}.}
    
    \item{version}{same as the \code{version} in 
    \code{\link{metaseqr2}}.}
    
    \item{db}{same as the \code{db} in 
    \code{\link{buildAnnotationDatabase}}.}
    
    \item{asdf}{return the result as a \code{\link{data.frame}}
    (default \code{FALSE}).}
    
    \item{rc}{same as the \code{rc} in 
    \code{\link{buildAnnotationDatabase}}.}
}
\value{
    The function returns a \code{GenomicRanges} object with
    the requested annotation.
}
\description{
    This function creates loads an annotation element from
    the local annotation database to be used with metaseqr2.
    If the annotation is not found and the organism is 
    supported, the annotation is created on the fly but not
    imported in the local database. Use
    \code{buildAnnotationDatabase} for this purpose.
    
}
\examples{
    db <- file.path(system.file(package="metaseqR2"),
        "annotation.sqlite")
    if (file.exists(db))
        gr <- loadAnnotation(genome="hg19", refdb="ensembl",
            level="gene", type="gene", db=db)
}
\author{
    Panagiotis Moulos
}

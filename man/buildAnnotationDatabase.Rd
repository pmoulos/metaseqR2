\name{buildAnnotationDatabase}
\alias{buildAnnotationDatabase}
\title{Build a local annotation database for metaseqR2}
\usage{
    buildAnnotationDatabase(organisms, sources,
    db = file.path(system.file(package = "metaseqR2"),
        "annotation.sqlite"),
    forceDownload = TRUE, rc = NULL)
}
\arguments{
    \item{organisms}{a list of organisms and versions
    for which to download and build annotations. Check 
    the main \code{\link{metaseqr2}} help page for details 
    on supported organisms and the Details section below.}

    \item{sources}{a character vector of public sources
    from which to download and build annotations. Check 
    the main \code{\link{metaseqr2}} help page for details 
    on supported annotation sources.}
    
    \item{db}{a valid path (accessible at least by the
    current user) where the annotation database will be 
    set up. It defaults to 
    \code{system.file(package = "metaseqR2"),}
    \code{"annotation.sqlite")} that is, the installation
    path of metaseqR2 package. See also Details.}
    
    \item{forceDownload}{by default, 
    \code{buildAnnotationDatabase} will not download an
    existing annotation again (\code{FALSE}). Set to 
    \code{TRUE} if you wish to update the annotation 
    database for a particular version.}

    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).
    Sometimes used for building certain annotation types.}
}
\value{
    The function does not return anything. Only the SQLite 
    database is created or updated.
}
\description{
    This function creates a local annotation database to be
    used with metaseqr2 so as to avoid long time on the fly 
    annotation downloads and formatting.
}
\details{
    Regarding the \code{organisms} argument, it is a list
    with specific format which instructs 
    \code{buildAnnotationDatabase} on which organisms and
    versions to download from the respective sources. Such
    a list may have the format: 
    \code{organisms=list(hg19=75, mm9=67, mm10=96:97)}
    This is explained as follows:
    \itemize{
        \item A database comprising the human genome versions
        \code{hg19} and the mouse genome versions 
        \code{mm9, mm10} will be constructed.
        \item If \code{"ensembl"} is in \code{sources}, 
        version 75 is downloaded for \code{hg19} and versions 
        \code{67, 96, 97} for \code{mm9, mm10}. 
        \item If \code{"ucsc"} or \code{"refseq"} are in 
        \code{sources}, the latest versions are downloaded
        and marked by the download date. As UCSC and RefSeq
        versions are not accessible in the same way as
        Ensembl, this procedure cannot always be replicated.
    }
    \code{organisms} can also be a character vector with organism
    names/versions (e.g. \code{organisms = c("mm10","hg19")}),
    then the latest versions are downloaded in the case of 
    Ensembl.
    
    Regarding \code{db}, this controls the location of the
    installation database. If the default is used, then there is
    no need to provide the local database path to any function
    that uses the database (e.g. the main \code{metaseqr2}).
    Otherwise, the user will either have to provide this each
    time, or the annotation will have to be downloaded and used
    on-the-fly.
}
\examples{
\donttest{
# Build a test database with one genome
myDb <- file.path(tempdir(),"testann.sqlite")

organisms <- list(mm9=67)
sources <- "ensembl"

# If the example is not running in a multicore system, rc is ignored
buildAnnotationDatabase(organisms,sources,db=myDb,rc=0.5)

# A more complete case, don't run as example
# Since we are using Ensembl, we can also ask for a version
#organisms <- list(
#    mm9=67,
#    mm10=96:97,
#    hg19=75,
#    hg38=96:97
#)
#sources <- c("ensembl", "refseq")

## Build on the default location (depending on package location, it may
## require root/sudo)
#buildAnnotationDatabase(organisms,sources)

## Build on an alternative location
#myDb <- file.path(path.expand("~"),"my_ann.sqlite")
#buildAnnotationDatabase(organisms,sources,db=myDb)
}
}
\author{
    Panagiotis Moulos
}

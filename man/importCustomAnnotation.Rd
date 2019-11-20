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
# Dummy GTF as example
chromInfo <- data.frame(length=c(1000L,2000L,1500L),
    row.names=c("A","B","C"))

# Build with the metadata list filled (you can also provide a version)
myGenes <- importCustomAnnotation(
  gtfFile=file.path(system.file(package="metaseqR2"),"dummy.gtf"),
  metadata=list(
    organism="dummy",
    source="dummy_db",
    version=1,
    chromInfo=chromInfo
  ),
  level="gene",type="gene"
)

## Real data!
## Gene annotation dump from Ensembl
#download.file(paste0("ftp://ftp.ensembl.org/pub/release-98/gtf/",
#  "dasypus_novemcinctus/Dasypus_novemcinctus.Dasnov3.0.98.gtf.gz"),
#  file.path(tempdir(),"Dasypus_novemcinctus.Dasnov3.0.98.gtf.gz"))

## Build with the metadata list filled (you can also provide a version)
#dasGenes <- importCustomAnnotation(
#  gtfFile=file.path(tempdir(),"Dasypus_novemcinctus.Dasnov3.0.98.gtf.gz"),
#  metadata=list(
#    organism="dasNov3_test",
#    source="ensembl_test"
#  ),
#  level="gene",type="gene"
#)
}
\author{
    Panagiotis Moulos
}

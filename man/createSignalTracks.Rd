\name{createSignalTracks}
\alias{createSignalTracks}
\title{Create bigWig signal tracks}
\usage{
    createSignalTracks(targets, org, urlBase = NULL, 
        stranded = FALSE, normTo = 1e+9, exportPath = ".",
        hubInfo = list(name = "MyHub", shortLabel = "My hub",
        longLabel = "My hub", email = "someone@example.com"),
        overwrite = FALSE, rc = NULL)
}
\arguments{
    \item{targets}{a tab-delimited file with the experimental 
    description or the output of \code{\link{readTargets}}.
    See also the \code{sampleList} argument in the main
    \code{\link{metaseqr2}} pipeline.}

    \item{org}{See the \code{org} argument in the main
    \code{\link{metaseqr2}} pipeline.}

    \item{urlBase}{a valid URL which is prepended to the created
    bigWig files.}

    \item{stranded}{Separate + and - strands and create separate
    bigWig files.}

    \item{normTo}{the total sum of signal to be used as the 
    normalization target. See also the \code{trackInfo} argument 
    in the main \code{\link{metaseqr2}} pipeline.}

    \item{exportPath}{path to export tracks.}
    
    \item{hubInfo}{information regarding the track hub created
    when \code{stranded=TRUE}. See also the \code{trackInfo} 
    argument in the main \code{\link{metaseqr2}} pipeline.}
    
    \item{overwrite}{overwrite tracks if they exist? Defaults to
    \code{FALSE}.}
    
    \item{rc}{Fraction of cores to use.}
}
\value{
    A string with the link(s) to the created tracks.
}
\description{
    This function creates bigWig files to be used for exploring
    RNA signal in genome browsers. When strands are separated,
    a UCSC genome browser trackhub is created to group tracks
    for the same sample. A link to the created data is returned.
}
\examples{
\dontrun{
# Will fail as BAM files do not exist
targets <- data.frame(samplename=c("C1","C2","T1","T2"),
  filename=c("C1_raw.bam","C2_raw.bam","T1_raw.bam","T2_raw.bam"),
  condition=c("Control","Control","Treatment","Treatment"),
  paired=rep("single",4),stranded=rep("forward",4))
path <- tempdir()
write.table(targets,file=file.path(path,"targets.txt"),
  sep="\t",row.names=FALSE,quote=FALSE)
link <- createSignalTracks(file.path(path,"targets.txt"),"hg19")
}
}
\author{
    Panagiotis Moulos
}

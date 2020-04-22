\name{read2count}
\alias{read2count}
\title{SAM/BAM/BED file reader helper for the metaseqr2 pipeline}
\usage{
    read2count(targets, annotation, fileType = targets$type,
        transLevel = "gene", utrOpts = list(frac = 1, 
        minLength = 300, downstream = 50), interFeature = FALSE, 
        rc = NULL)
}
\arguments{
    \item{targets}{a named list, the output of 
    \code{\link{readTargets}} or an existing file with
    targets. See also the main \code{metaseqr2} man page.}

    \item{annotation}{a \code{GenomicRanges} or 
    \code{data.frame} with genomic coordinates to use for
    read counting. See also \code{\link{getAnnotation}}.}

    \item{fileType}{the type of raw input files. It can be
    \code{"bed"} for BED files or \code{"sam"}, \code{"bam"}
    for SAM/BAM files. See the same argument in the main
    \code{\link{metaseqr2}} function for the case of
    auto-guessing.}
    
    \item{transLevel}{see the \code{transLevel} argument
    in the main \code{\link{metaseqr2}} function.}
    
    \item{utrOpts}{a named list with members \code{frac} 
    which is the fraction (0-1) of the 3' UTR region to count 
    reads in, \code{minLength} the minimum acceptable 3'UTR
    length irrespective of \code{frac} and \code{downstream} 
    the number of base pairs to flank the end of the 3' UTR of 
    transcripts when analyzing Quant-Seq data.}
    
    \item{interFeature}{see the \code{inter.feature} argument
    in \code{summarizeOverlaps}.}

    \item{rc}{the fraction of the available cores to use
    in a multicore system.}
}
\value{
    A data frame with counts for each sample, ready to be
    passed to the main \code{\link{metaseqr2}} pipeline.
}
\description{
    This function is a helper for the \code{metaseqr2}
    pipeline, for reading SAM/BAM or BED files when a read
    counts file is not available. It can also be used
    very easily in an autonomous manner.
}
\examples{
dataPath <- system.file("extdata",package="metaseqR2")
targets <- data.frame(samplename=c("C","T"),
  filename=file.path(dataPath,c("C.bam","T.bam")),  
  condition=c("Control","Treatment"),
  paired=c("single","single"),stranded=c("forward","forward"))
path <- tempdir()
write.table(targets,file=file.path(path,"targets.txt"),
  sep="\t",row.names=FALSE,quote=FALSE)
geneData <- loadAnnotation("mm9","ensembl","gene")
myTargets <- readTargets(file.path(path,"targets.txt"))
if (.Platform$OS.type == "unix") {
  r2c <- read2count(targets=myTargets,
    fileType=myTargets$type,annotation=geneData)
  geneCounts <- r2c$counts
  libsizeList <- r2c$libsize
}
}
\author{
    Panagiotis Moulos
}


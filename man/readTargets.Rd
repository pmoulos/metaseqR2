\name{readTargets}
\alias{readTargets}
\title{Creates sample list and BAM/BED file list from file}
\usage{
    readTargets(input, path = NULL)
}
\arguments{
    \item{input}{a tab-delimited file specifically
    structured. See Details.}

    \item{path}{an optional path where all the BED/BAM 
    files are placed, to be prepended to the BAM/BED 
    file names in the targets file.}
}
\value{
    A named list with four members. The first member is 
    a named list whose names are the conditions of the 
    experiments and its members are the samples belonging 
    to each condition. The second member is like the 
    first, but this time the members are named vectors 
    whose names are the sample names and the vector 
    elements are full path to BAM/BED files. The third 
    member is like the second, but instead of filenames 
    it contains information about single- or paired-end
    reads (if available). The fourth member is like the 
    second, but instead of filenames it contains 
    information about the strandedness of the reads (if
    available). The fifth member is the guessed type 
    of the input files (SAM/BAM or BED). It will be used 
    if not given in the main \code{\link{read2count}}
    function.
}
\description{
    Create the main sample list and determine the BAM/BED
    files for each sample from an external file.
}
\details{
	Regarding the input file, its files must be structured 
    as follows: the first line of the external tab delimited
    file should contain column names (names are not
    important). The first column MUST contain UNIQUE sample 
    names. The second column MUST contain the raw BAM/BED 
    files WITH their full path. Alternatively, the \code{path} 
    argument should be provided (see below). The third column 
    MUST contain the biological condition where each of the 
    samples in the first column should belong to. There is an 
    optional fourth column which should contain the keywords 
    \code{"single"} for single-end reads, \code{"paired"} for 
    paired-end reads or \code{"mixed"} for BAM files that 
    contain both single- and paired-end reads (e.g. after a 
    mapping procedure with two round of alignment). If this 
    column is not provided, single-end reads will be assumed. 
    There is an optional fifth column which stranded read 
    assignment. It should contain the keywords \code{"forward"} 
    for a forward (5'->3') strand library construction protocol,
    \code{"reverse"} for a reverse (3'->5') strand library 
    construction protocol, or \code{"no"} for unstranded/unknown 
    protocol. If this column is not provided, unstranded reads 
    will be assumed.
}
\examples{
targets <- data.frame(sample=c("C1","C2","T1","T2"),
    filename=c("C1_raw.bam","C2_raw.bam","T1_raw.bam","T2_raw.bam"),
    condition=c("Control","Control","Treatment","Treatment"))
path <- tempdir()
write.table(targets,file=file.path(path,"targets.txt"),
    sep="\\t",row.names=FALSE,quote=FALSE)
theList <- readTargets(file.path(path,"targets.txt"),path=path)
sampleList <- theList$samples
bamfileList <- theList$files
}
\author{
    Panagiotis Moulos
}

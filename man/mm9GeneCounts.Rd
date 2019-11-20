\docType{data}
\name{mm9GeneCounts}
\alias{mm9GeneCounts}
\title{mouse RNA-Seq data with two conditions, four samples}
\format{a \code{data.frame} with gene read counts and some embedded 
annotation, one row per gene.}
\source{
    ENCODE (http://genome.ucsc.edu/encode/)
}
\description{
    This data set contains RNA-Seq gene read counts for 3
    chromosomes. The data were downloaded from the ENCODE
    public repository and are derived from the study of
    Mortazavi et al., 2008 (Mortazavi A, Williams BA, McCue
    K, Schaeffer L, Wold B. Mapping and quantifying mammalian
    transcriptomes by RNA-Seq. Nat Methods. 2008
    Jul;5(7):621-8). In their experiment, the authors studied
    among others genes expression at two developmental stages
    of mouse liver cells. It has two conditions-developmental
    stages (e14.5, adult_8_weeks) and four samples (e14.5_1,
    e14.5_2, a8w_1, a8w_2). It also contains a predefined
    \code{sampleList} and \code{libsizeList} named
    \code{sampleListMm9} and \code{libsizeListMm9}.
}
\author{
    Panagiotis Moulos
}
\keyword{datasets}

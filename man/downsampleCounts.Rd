\name{downsampleCounts}
\alias{downsampleCounts}
\title{Downsample read counts}
\usage{
    downsampleCounts(counts,seed=42)
}
\arguments{
    \item{counts}{the read counts table 
    which is subjected to downsampling.}

    \item{seed}{random seed for reproducible 
    downsampling.}
}
\value{
    The downsampled counts matrix.
}
\description{
    This function downsamples the library sizes 
    of a read counts table to the lowest library 
    size, according to the methdology used in 
    (Soneson and Delorenzi, BMC Bioinformatics, 
    2013).
}
\examples{
require(DESeq)
dataMatrix <- counts(makeExampleCountDataSet())
D <- downsampleCounts(dataMatrix)
}
\author{
    Panagiotis Moulos
}

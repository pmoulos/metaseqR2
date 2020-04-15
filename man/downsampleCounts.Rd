\name{downsampleCounts}
\alias{downsampleCounts}
\title{Downsample read counts}
\usage{
    downsampleCounts(counts)
}
\arguments{
    \item{counts}{the read counts table 
    which is subjected to downsampling.}
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
\details{
    The downsampling process involves random sampling. 
    For guaranteed reproducibility, be sure to use 
    \code{set.seed} before downsampling. By default, 
    when the metaseqR2 package is loaded, the seed is 
    set to \code{42}.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(5000)
D <- downsampleCounts(dataMatrix)
}
\author{
    Panagiotis Moulos
}

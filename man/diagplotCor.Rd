\name{diagplotCor}
\alias{diagplotCor}
\title{Summarized correlation plots}
\usage{
    diagplotCor(mat, type = c("heatmap", "correlogram"),
        output = "x11", path = NULL, ...)
}
\arguments{
    \item{mat}{the read counts matrix or data frame.}

    \item{type}{create heatmap of correlogram plots.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filename of the pairwise comparisons plot produced if
    it's a file.
}
\description{
    This function uses the read counts matrix to create
    heatmap or correlogram correlation plots.
}
\examples{
require(DESeq)
dataMatrix <- counts(makeExampleCountDataSet())
diagplotCor(dataMatrix,type="heatmap")
diagplotCor(dataMatrix,type="correlogram")
}
\author{
    Panagiotis Moulos
}


\name{diagplotMds}
\alias{diagplotMds}
\title{Multi-Dimensinal Scale plots or RNA-Seq samples}
\usage{
    diagplotMds(x, sampleList, method = "spearman",
        logIt = TRUE, output = "x11", path = NULL, ...)
}
\arguments{
    \item{x}{the count data matrix.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{method}{which correlation method to use. Same as
    the method parameter in \code{\link{cor}} function.}

    \item{logIt}{whether to log transform the values of x or
    not.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"}, \code{"ps"} or \code{"json"}. The latter is
    currently available for the creation of interactive
    volcano plots only when reporting the output, through the
    highcharts javascript library.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filename of the MDS plot produced if it's a file.
}
\description{
    Creates a Multi-Dimensional Scale plot for the given
    samples based on the count data matrix. MDS plots are
    very useful for quality control as you can easily see of
    samples of the same groups are clustered together based
    on the whole dataset.
}
\examples{
require(DESeq)
dataMatrix <- counts(makeExampleCountDataSet())
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplotMds(dataMatrix,sampleList)
}
\author{
    Panagiotis Moulos
}


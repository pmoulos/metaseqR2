\name{diagplotFiltered}
\alias{diagplotFiltered}
\title{Diagnostic plot for filtered genes}
\usage{
    diagplotFiltered(x, y, output = "x11", path = NULL, ...)
}
\arguments{
    \item{x}{an annotation data frame like the ones produced
    by \code{\link{get.annotation}}. \code{x} should be the
    filtered annotation according to metaseqR's filters.}

    \item{y}{an annotation data frame like the ones produced
    by \code{\link{get.annotation}}. \code{y} should contain
    the total annotation without the application of any
    metaseqr filter.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filenames of the plots produced in a named list with
    names the \code{which.plot} argument. If
    output=\code{"x11"}, no output filenames are produced.
}
\description{
    This function plots a grid of four graphs depicting: in
    the first row, the numbers of filtered genes per
    chromosome in the first column and per biotype in the
    second column. In the second row, the percentages of
    filtered genes per chromosome related to the whole genome
    in the first columns and per biotype in the second
    column.
}
\examples{
data("mm9GeneData",package="metaseqR2")
y <- mm9GeneCounts[,c(1:6,8,7)]
x <- y[-sample(1:nrow(y),1000),]
diagplotFiltered(x,y)
}
\author{
    Panagiotis Moulos
}

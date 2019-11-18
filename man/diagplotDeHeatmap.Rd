\name{diagplotDeHeatmap}
\alias{diagplotDeHeatmap}
\title{Diagnostic heatmap of differentially expressed genes}
\usage{
    diagplotDeHeatmap(x, scale = c("asis", "zscore"), con = NULL, 
        output = "x11", path = NULL, ...)
}
\arguments{
    \item{x}{the data matrix to create a heatmap for.}
    
    \item{scale}{value scale in the heatmap. As provided 
    (\code{scale="asis"}, default) or Z-scores
    (\code{scale="zscore"})}

    \item{con}{an optional string depicting a name (e.g. the
    contrast name) to appear in the title of the volcano
    plot.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"}, \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filenames of the plots produced in a named list with
    names the \code{whichPlot} argument. If
    \code{output="x11"}, no output filenames are produced.
}
\description{
    This function plots a heatmap of the differentially
    expressed genes produced by the metaseqr workflow, useful
    for quality control, e.g. whether samples belonging to
    the same group cluster together.
}
\examples{
require(DESeq)
dataMatrix <- counts(makeExampleCountDataSet())
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
M <- normalizeEdger(dataMatrix,sampleList)
p <- statEdger(M,sampleList,contrast)
diagplotDeHeatmap(dataMatrix[p[[1]]<0.05,])
}
\author{
    Panagiotis Moulos
}


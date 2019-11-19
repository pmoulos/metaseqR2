\name{diagplotNoiseq}
\alias{diagplotNoiseq}
\title{Diagnostic plots based on the NOISeq package}
\usage{
    diagplotNoiseq(x, sampleList, covars,
        whichPlot = c("biodetection", "countsbio", "saturation", 
        "rnacomp", "readnoise", "biodist"),
        output = "x11",
        biodistOpts = list(p = NULL, pcut = NULL, name = NULL),
        path = NULL, isNorm = FALSE, ...)
}
\arguments{
    \item{x}{the count data matrix.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{covars}{a list (whose annotation elements are
    ideally a subset of an annotation data frame produced by
    \code{\link{getAnnotation}}) with the following members:
    data (the data matrix), length (gene length), gc (the
    gene gc_content), chromosome (a data frame with
    chromosome name and co-ordinates), factors (a factor with
    the experimental condition names replicated by the number
    of samples in each experimental condition) and biotype
    (each gene's biotype as depicted in Ensembl-like
    annotations).}

    \item{whichPlot}{the NOISeq package plot to generate. 
    See Details}

    \item{biodistOpts}{a list with the following members: p
    (a vector of p-values, e.g. the p-values of a contrast),
    pcut (a unique number depicting a p-value cutoff,
    required for the \code{"biodist"} case), name (a name for
    the \code{"biodist"} plot, e.g. the name of the
    contrast.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{isNorm}{a logical indicating whether object
    contains raw or normalized data. It is not essential and
    it serves only plot annotation purposes.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filenames of the plots produced in a named list with
    names the \code{whichPlot} argument. If
    \code{output="x11"}, no output filenames are produced.
}
\description{
    A wrapper around the plotting functions availale in the
    NOISeq Bioconductor package. For analytical explanation
    of each plot please see the vignette of the NOISeq 
    package. It is best to use this function through the 
    main plotting function \code{\link{metaseqrPlot}}.
}
\details{
	Regarding \code{whichPlot}, It can be one or more of 
	\code{"biodetection"}, \code{"countsbio"}, 
	\code{"saturation"}, \code{"rnacomp"}, \code{"readnoise"} 
	or \code{"biodist"}. Please refer to the documentation of 
	the NOISeq package for details on the use of these plots. 
	The \code{whichPlot="saturation"} case is modified to be
    more informative by producing two kinds of plots.
}
\note{
    Please note that in case of \code{"biodist"} plots, the
    behavior of the function is unstable, mostly due to the
    very specific inputs this plotting function accepts in
    the NOISeq package. We have tried to predict unstable
    behavior and avoid exceptions through the use of tryCatch
    but it's still possible that you might run onto an error.
}
\examples{
require(DESeq)
dataMatrix <- counts(makeExampleCountDataSet())
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
lengths <- round(1000*runif(nrow(dataMatrix)))
starts <- round(1000*runif(nrow(dataMatrix)))
ends <- starts + lengths
covars <- list(
  data=dataMatrix,
  length=lengths,
  gc=runif(nrow(dataMatrix)),
  chromosome=data.frame(
    chromosome=c(rep("chr1",nrow(dataMatrix)/2),
      rep("chr2",nrow(dataMatrix)/2)),
    start=starts,
    end=ends
  ),
  factors=data.frame(class=metaseqR2:::asClassVector(sampleList)),
  biotype=c(rep("protein_coding",nrow(dataMatrix)/2),rep("ncRNA",
    nrow(dataMatrix)/2))
)
p <- runif(nrow(dataMatrix))
diagplotNoiseq(dataMatrix,sampleList,covars=covars,
  biodistOpts=list(p=p,pcut=0.1,name="A_vs_B"))
}
\author{
    Panagiotis Moulos
}

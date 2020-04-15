\name{metaseqrPlot}
\alias{metaseqrPlot}
\title{Diagnostic plots for the metaseqR2 package}
\usage{
    metaseqrPlot(object, sampleList, annotation = NULL,
        contrastList = NULL, pList = NULL,
        thresholds = list(p = 0.05, f = 1),
        plotType = c("mds", "biodetection", "countsbio",
            "saturation", "readnoise", "rnacomp", "correl",
            "pairs", "boxplot", "gcbias", "lengthbias",
            "meandiff", "meanvar", "deheatmap", "volcano",
            "biodist", "filtered", "mastat", "deregulogram",
            "statvenn", "foldvenn"),
        isNorm = FALSE, output = "x11", path = NULL, ...)
}
\arguments{
    \item{object}{a matrix or a data frame containing count
    data derived before or after the normalization procedure,
    filtered or not by the metaseqR2's filters and/or p-value.}

    \item{sampleList}{the list containing condition names
    and the samples under each condition.}

    \item{annotation}{a data frame containing annotation
    elements for each row in object. See also Details.}

    \item{contrastList}{the vector of contrasts as defined in 
    the main help page of \code{\link{metaseqr2}}.}

    \item{pList}{a list of p-values for each contrast as
    obtained from any of the \code{stat*} methods of the
    metaseqr package. See also Details.}

    \item{thresholds}{a list with the elements \code{"p"} and
    \code{"f"} which are the p-value and the fold change
    cutoff when \code{diagplotType="volcano"}.}

    \item{plotType}{one or more of the diagnostic plots
    supported in metaseqR2 package. See also Details.}

    \item{isNorm}{a logical indicating whether object
    contains raw or normalized data. It is not essential and
    it serves only plot annotation purposes.}

    \item{output}{one or more R plotting device to direct the
    plot result to. See Details.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    A named list containing the file names of the produced
    plots. Each list member is names according to the
    selected plotting device and is also a named list, whose
    names are the plot types. The final contents are the file
    names in case the plots are written to a physical
    location (not meaningful for \code{"x11"}).
}
\description{
    This is the main function for producing sructured quality
    control and informative graphs base on the results of the
    various steps of the metaseqR package. The graphs
    produced span a variety of issues like good sample
    reproducibility (Multi-Dimensional Scaling plot, biotype
    detection, heatmaps. diagplotMetaseqr, apart from
    implementing certain package-specific plots, is a wrapper
    around several diagnostic plots present in other RNA-Seq
    analysis packages such as EDASeq and NOISeq.
}
\details{
    Regarding \code{object}, the object can be fed to any of 
    the \code{diagplotMetaseqr} plotting systems but not every
    plot is meaningful. For example, it's meaningless to
    create a \code{"biodist"} plot for a count matrix before
    normalization or statistical testing.
    
    Regarding \code{annotation}, usually, it is a subset of 
    the annotation obtained by \code{\link{getAnnotation}} or
    a subset of possibly embedded annotation with the input
    counts table. This parameter is optional and required
    only when diagplotType is any of \code{"biodetection"},
    \code{"countsbio"}, \code{"saturation"},
    \code{"rnacomp"}, \code{"readnoise"}, \code{"biodist"},
    \code{"gcbias"}, \code{"lengthbias"} or
    \code{"filtered"}.
    
    Regarding \code{contrastList}, this parameter is optional 
    and required only when \code{diagplotType} is any of
    \code{"deheatmap"}, \code{"volcano"} or \code{"biodist"}.
    It can also be a named structured list of contrasts as 
    returned by the internal function 
    \code{metaseqR2:::makeContrastList}.
    
    Regarding \code{diagplotType}, many of these plots
    require the presence of additional package, something
    that is checked while running the main metaseqr2 function.
    The supported plots are \code{"mds"}, \code{"biodetection"},
    \code{"countsbio"}, \code{"saturation"}, \code{"rnacomp"},
    \code{"boxplot"}, \code{"gcbias"}, \code{"lengthbias"},
    \code{"meandiff"}, \code{"meanvar"}, \code{"deheatmap"},
    \code{"volcano"}, \code{"biodist"}, \code{"filtered"}, 
    \code{"readnoise"}, \code{"venn"}, \code{"correl"},
    \code{"pairwise"}. For a brief description of these plots 
    please see the main \code{\link{metaseqr2}} help page.
    
    Regarding \code{pList}, this parameter is optional 
    and required only when \code{diagplotType} is any of
    \code{"deheatmap"}, \code{"volcano"} or \code{"biodist"}.
    
    Regarding \code{output}, supported mechanisms are: 
    \code{"png"}, \code{"jpg"}, \code{"bmp"}, \code{"pdf"},
    \code{"ps"} or \code{"json"}. The latter is currently 
    available for the creation of interactive volcano plots 
    only when reporting the output, through the highcharts 
    javascript library. The default plotting (\code{"x11"}) 
    is not supported due to instability in certain devices.
}
\note{
    In order to make the best out of this function, you
    should generally provide the annotation argument as most
    and also the most informative plots depend on this. If
    you don't know what is inside your counts table or how
    many annotation elements you can provide by embedding it,
    it's always best to setup a local databse so as to use
    predefined annotations that work better with the
    functions of the whole package.
}
\examples{
dataMatrix <- metaseqR2:::exampleCountData(2000)
sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
metaseqrPlot(dataMatrix,sampleList,plotType=c("mds","boxplot"))

normArgs <- getDefaults("normalization","deseq2")
object <- normalizeDeseq2(dataMatrix,sampleList,normArgs)
metaseqrPlot(object,sampleList,plotType="boxplot")

## More
#p <- statDeseq2(object,sampleList)
#metaseqrPlot(object,sampleList,contrastList=contrast,pList=p,
#    plotType="volcano")
}
\author{
    Panagiotis Moulos
}

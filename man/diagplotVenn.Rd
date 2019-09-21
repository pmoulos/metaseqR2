\name{diagplotVenn}
\alias{diagplotVenn}
\title{Venn diagrams when performing meta-analysis}
\usage{
    diagplotVenn(pmat, fcmat = NULL, pcut = 0.05,
        fcut = 0.5, direction = c("dereg", "up", "down"),
        nam = as.character(round(1000 * runif(1))),
        output = "x11", path = NULL, altNames = NULL, ...)
}
\arguments{
    \item{pmat}{a matrix with p-values corresponding to the
    application of each statistical algorithm. See also
    Details.}

    \item{fcmat}{an optional matrix with fold changes
    corresponding to the application of each statistical
    algorithm. See also Details.}

    \item{pcut}{if \code{fcmat} is supplied, an absolute
    fold change cutoff to be applied to \code{fcmat} to
    determine the differentially expressed genes for each
    algorithm.}

    \item{fcut}{a p-value cutoff for statistical
    significance. Defaults to \code{0.05}.}

    \item{direction}{if \code{fcmat} is supplied, a keyword
    to denote which genes to draw in the Venn diagrams with
    respect to their direction of regulation. See Details.}

    \item{nam}{a name to be appended to the output graphics
    file (if \code{"output"} is not \code{"x11"}).}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files. If 
    \code{"path"} is not \code{NULL}, a file with the 
    intersections in the Venn diagrams will be produced 
    and written in \code{"path"}.}

    \item{altNames}{an optional named vector of names, e.g.
    HUGO gene symbols, alternative or complementary to the
    unique gene names which are the rownames of \code{pmat}.
    The names of the vector must be the rownames of
    \code{pmat}.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filenames of the plots produced in a named list with
    names the \code{which.plot} argument. If
    output=\code{"x11"}, no output filenames are produced.
}
\description{
    This function uses the R package VennDiagram and plots an
    up to 5-way Venn diagram depicting the common and
    specific to each statistical algorithm genes, for each
    contrast. Mostly for internal use because of its main
    argument which is difficult to construct, but can be used
    independently if the user grasps the logic.
}
\details{
	Regarding \code{pmat}, the p-value matrix must have the 
	colnames attribute and the colnames should correspond to 
	the name of the algorithm used to fill the specific 
	column (e.g. if \code{"statistics"=c("deseq","edger",}
	\code{"nbpseq")} then \code{colnames(pmat) <-}
    \code{c("deseq","edger","nbpseq")}.
    
    Regarding \code{fcmat}, the fold change matrix must have 
    the colnames attribute and the colnames should correspond 
    to the name of the algorithm used to fill the specific 
    column (see the parameter \code{pmat}).
    
    Regarding \code{direction}, it can be one of \code{"dereg"}
    for the total of regulated genes, where 
    \code{abs(fcmat[,n])>=fcut} (default), \code{"up"} for
    the up-regulated genes where \code{fcmat[,n]>=fcut} or
    \code{"down"} for the up-regulated genes where
    \code{fcmat[,n]<=-fcut}.
}
\examples{
require(VennDiagram)
p1 <- 0.01*matrix(runif(300),100,3)
p2 <- matrix(runif(300),100,3)
p <- rbind(p1,p2)
rownames(p) <- paste("gene",1:200,sep="_")
colnames(p) <- paste("method",1:3,sep="_")
vennContents <- diagplotVenn(p)
}
\author{
    Panagiotis Moulos
}

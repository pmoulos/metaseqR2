\name{diagplotFtd}
\alias{diagplotFtd}
\title{Create False (or True) Positive (or
    Negative) curves}
\usage{
    diagplotFtd(truth, p, type = "fpc", N = 2000,
        output = "x11", path = NULL, draw = TRUE, ...)
}
\arguments{
    \item{truth}{the ground truth differential 
    expression vector. It should contain only 
    zero and non-zero elements, with zero denoting
    non-differentially expressed genes and non-zero, 
    differentially expressed genes. Such a vector 
    can be obtained for example by using the 
    \code{\link{makeSimDataSd}} function, which 
    creates simulated RNA-Seq read counts based on 
    real data. The elements of \code{truth} MUST 
    be named (e.g. each gene's name).}

    \item{p}{a p-value matrix whose rows correspond 
    to each element in the \code{truth} vector. If 
    the matrix has a \code{colnames} attribute, a 
    legend will be added to the plot using these 
    names, else a set of column names will be 
    auto-generated. \code{p} can also be a list or 
    a data frame. The p-values MUST be named (e.g. 
    each gene's name).}

    \item{type}{what to plot, can be \code{"fpc"} 
    for False Positive Curves (default), 
    \code{"tpc"} for True Positive Curves, 
    \code{"fnc"} for False Negative Curves or 
    \code{"tnc"} for True Negative Curves.}

    \item{N}{create the curves based on the 
    top (or bottom) \code{N} ranked genes 
    (default is 2000) to be used with
    \code{type="fpc"} or \code{type="tpc"}.}

    \item{output}{one or more R plotting device to 
    direct the plot result to. Supported mechanisms: 
    \code{"x11"} (default), \code{"png"}, \code{"jpg"}, 
    \code{"bmp"}, \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{draw}{boolean to determine whether to
    plot the curves or just return the calculated
    values (in cases where the user wants the
    output for later averaging for example). 
    Defaults to \code{TRUE} (make plots).}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    A named list with two members: the first member 
    (\code{ftdr}) contains the values used to create 
    the plot. The second member (\code{path}) contains 
    the path to the created figure graphic.
}
\description{
    This function creates false (or true) discovery 
    curves using a matrix of p-values (such a matrix 
    can be derived for example from the result table 
    of \code{\link{metaseqr2}} by subsetting the table 
    to get the p-values from several algorithms) 
    given a ground truth vector for differential 
    expression.
}
\examples{
p1 <- 0.001*matrix(runif(300),100,3)
p2 <- matrix(runif(300),100,3)
p <- rbind(p1,p2)
rownames(p) <- paste("gene",1:200,sep="_")
colnames(p) <- paste("method",1:3,sep="_")
truth <- c(rep(1,40),rep(-1,40),rep(0,20),
    rep(1,10),rep(2,10),rep(0,80))
names(truth) <- rownames(p)
ftdObj <- diagplotFtd(truth,p,N=100)
}
\author{
    Panagiotis Moulos
}


\name{diagplotRoc}
\alias{diagplotRoc}
\title{Create basic ROC curves}
\usage{
    diagplotRoc(truth, p, sig = 0.05, x = "fpr", 
        y = "tpr", output = "x11", path = NULL,
        draw = TRUE, ...)
}
\arguments{
    \item{truth}{the ground truth differential 
    expression vector. It should contain only 
    zero and non-zero elements, with zero denoting
    non-differentially expressed genes and non-zero, 
    differentially expressed genes. Such a vector 
    can be obtained for example by using the 
    \code{\link{make.sim.data.sd}} function, which 
    creates simulated RNA-Seq read counts based on 
    real data.}

    \item{p}{a p-value matrix whose rows correspond 
    to each element in the \code{truth} vector. If 
    the matrix has a \code{colnames} attribute, a 
    legend will be added to the plot using these 
    names, else a set of column names will be 
    auto-generated. \code{p} can also be a list or 
    a data frame.}

    \item{sig}{a significance level (0 < \code{sig} 
    <=1).}

    \item{x}{what to plot on x-axis, can be one of 
    \code{"fpr"}, \code{"fnr"}, \code{"tpr"}, 
    \code{"tnr"} for False Positive Rate, False
    Negative Rate, True Positive Rate and True 
    Negative Rate respectively.}

    \item{y}{what to plot on y-axis, same as 
    \code{x} above.}

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
    A named list with two members. The first member is 
    a list containing the ROC statistics: \code{TP} 
    (True Postives), \code{FP} (False Positives), 
    \code{FN} (False Negatives), \code{TN} 
    (True Negatives), \code{FPR} (False Positive Rate), 
    \code{FNR} (False Negative Rate), \code{TPR} (True 
    Positive Rate), \code{TNR} (True Negative Rate), 
    \code{AUC} (Area Under the Curve). The second is 
    the path to the created figure graphic.
}
\description{
    This function creates basic ROC curves using a 
    matrix of p-values (such a matrix can be 
    derived for example from the result table of 
    \code{\link{metaseqr}} by subsetting the table 
    to get the p-values from several algorithms) 
    given a ground truth vector for differential 
    expression and a significance level.
}
\examples{
p1 <- 0.001*matrix(runif(300),100,3)
p2 <- matrix(runif(300),100,3)
p <- rbind(p1,p2)
rownames(p) <- paste("gene",1:200,sep="_")
colnames(p) <- paste("method",1:3,sep="_")
truth <- c(rep(1,40),rep(-1,40),rep(0,20),rep(1,10),
    rep(2,10),rep(0,80))
names(truth) <- rownames(p)
rocObj <- diagplotRoc(truth,p)
}
\author{
    Panagiotis Moulos
}


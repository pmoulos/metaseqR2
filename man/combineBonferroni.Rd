\name{combineBonferroni}
\alias{combineBonferroni}
\title{Combine p-values with Bonferroni's method}
\usage{
    combineBonferroni(p, zerofix = NULL)
}
\arguments{
    \item{p}{a vector of p-values for each statistical 
    tests).}
    
    \item{zerofix}{\code{NULL} (default) or a fixed 
    numeric value between 0 and 1.}
}
\value{
    A vector of combined p-values. 
}
\description{
    This function combines p-values from the 
    various statistical tests supported by 
    metaseqR2 using the Bonferroni's method (see 
    reference in the main \code{\link{metaseqr2}}
    help page or in the vignette).
}
\details{
    The argument \code{zerofix} is used to correct for
    the case of a p-value which is equal to 0 as a result
    of internal numerical and approximation procedures.
    When \code{NULL}, random numbers greater than 0 and
    less than or equal to 0.5 are used to multiply the
    offending p-values with the lowest provided non-zero
    p-value, maintaining thus a virtual order of 
    significance, avoiding having the same p-values for 
    two tests and assuming that all zero p-values represent
    extreme statistical significance. When a numeric
    between 0 and 1, this number is used for the above
    multiplication instead.
}
\examples{
p <- matrix(runif(300),100,3)
pc <- combineBonferroni(p)
}
\author{
    Panagiotis Moulos
}


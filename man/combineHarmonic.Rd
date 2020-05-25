\name{combineHarmonic}
\alias{combineHarmonic}
\title{Combine p-values using weights}
\usage{
    combineHarmonic(p, w, zerofix = NULL)
}
\arguments{
    \item{p}{a p-value matrix (rows are genes, 
    columns are statistical tests).}

    \item{w}{a weights vector, must sum to 1.}
    
    \item{zerofix}{\code{NULL} (default) or a fixed 
    numeric value between 0 and 1.}
}
\value{
    A vector of combined p-values. 
}
\description{
    This function combines p-values from the 
    various statistical tests supported by
    metaseqR using p-value weights.
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
pc <- combineHarmonic(p,w=c(0.2,0.5,0.3))
}
\author{
    Panagiotis Moulos
}

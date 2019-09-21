\name{combineWeight}
\alias{combineWeight}
\title{Combine p-values using weights}
\usage{
    combineWeight(p, w)
}
\arguments{
    \item{p}{a p-value matrix (rows are genes, 
    columns are statistical tests).}

    \item{w}{a weights vector, must sum to 1.}
}
\value{
    A vector of combined p-values. 
}
\description{
    This function combines p-values from the 
    various statistical tests supported by
    metaseqR using p-value weights.
}
\examples{
p <- matrix(runif(300),100,3)
pc <- combineWeight(p,w=c(0.2,0.5,0.3))
}
\author{
    Panagiotis Moulos
}


\name{combineSimes}
\alias{combineSimes}
\title{Combine p-values with Simes' method}
\usage{
    combineSimes(p)
}
\arguments{
    \item{p}{a p-value matrix (rows are genes, 
    columns are statistical tests).}
}
\value{
    A vector of combined p-values. 
}
\description{
    This function combines p-values from the 
    various statistical tests supported by 
    metaseqR using the Simes' method (see 
    reference in the main \code{\link{metaseqr2}}
    help page or in the vignette).
}
\examples{
p <- matrix(runif(300),100,3)
pc <- combineSimes(p)
}
\author{
    Panagiotis Moulos
}


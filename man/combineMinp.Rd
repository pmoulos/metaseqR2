\name{combineMinp}
\alias{combineMinp}
\title{Combine p-values using the minimum p-value}
\usage{
    combineMinp(p)
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
    metaseqR by taking the minimum p-value.
}
\examples{
p <- matrix(runif(300),100,3)
pc <- combineMinp(p)
}
\author{
    Panagiotis Moulos
}


\name{combineMaxp}
\alias{combineMaxp}
\title{Combine p-values using the maximum p-value}
\usage{
    combineMaxp(p)
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
    metaseqR by taking the maximum p-value.
}
\examples{
p <- matrix(runif(300),100,3)
pc <- combineMaxp(p)
}
\author{
    Panagiotis Moulos
}


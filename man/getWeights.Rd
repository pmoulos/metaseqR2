\name{getWeights}
\alias{getWeights}
\title{Get precalculated statistical test 
    weights}
\usage{
    getWeights(org = c("human", "chimpanzee", "mouse",
        "fruitfly", "arabidopsis", "rat"))
}
\arguments{
    \item{org}{\code{"human"}, \code{"chimpanzee"},
    \code{"mouse"}, \code{"fruitfly"}, \code{"arabidopsis"}
    or \code{"rat"}.}
}
\value{
    A named vector of convex weights.
}
\description{
    This function returns pre-calculated weights 
    for human, chimpanzee, mouse, fruitfly and 
    arabidopsis based on the performance of 
    simulated datasets estimated from real data 
    from the ReCount database 
    (\url{http://bowtie-bio.sourceforge.net/recount/}).
    Currently pre-calculated weights are available only 
    when all six statistical tests are used and for 
    normalization with EDASeq. For other combinations, 
    use the \code{\link{estimateAufcWeights}} function.
}
\examples{
wh <- getWeights("human")
}
\author{
    Panagiotis Moulos
}

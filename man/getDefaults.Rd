\name{getDefaults}
\alias{getDefaults}
\title{Default parameters for several metaseqr functions}
\usage{
    getDefaults(what, method = NULL)
}
\arguments{
    \item{what}{a keyword determining the procedure for which
    to fetch the default settings according to method
    parameter. It can be one of \code{"normalization"},
    \code{"statistics"}, \code{"geneFilter"},
    \code{"exonFilter"} or \code{"biotypeFilter"}.}

    \item{method}{the supported algorithm included in
    metaseqR for which to fetch the default settings.
    Se Details.}
}
\value{
    A list with default setting that can be used directly in
    the call of metaseqr.
}
\description{
    This function returns a list with the default settings
    for each filtering, statistical and normalization
    algorithm included in the metaseqR package. See the
    documentation of the main function and the documentation
    of each statistical and normalization method for details.
}
\details{
    When \code{what} is \code{"normalization"}, method is 
    one of \code{"edaseq"}, \code{"deseq"}, \code{"edger"},
    \code{"noiseq"} or \code{"nbpseq"}. When \code{what} is
    \code{"statistics"}, method is one of \code{"deseq"},
    \code{"edger"}, \code{"noiseq"}, \code{"bayseq"},
    \code{"limma"} or \code{"nbpseq"}. When \code{method} is
    \code{"biotypeFilter"}, \code{what} is the input
    organism (see the main \code{\link{metaseqr2}} help page
    for a list of supported organisms).
}
\examples{
normArgsEdaseq <- getDefaults("normalization","edaseq")
statArgsEdger <- getDefaults("statistics","edger")
}
\author{
    Panagiotis Moulos
}

\name{getInstalledAnnotations}
\alias{getInstalledAnnotations}
\title{Load a metaseqR2 annotation element}
\usage{
    getInstalledAnnotations(obj = NULL)
}
\arguments{
    \item{obj}{\code{NULL} or the path to a metaseqR2
    SQLite annotation database. If \code{NULL}, the
    function will try to guess the location of the 
    SQLite database.}
}
\value{
    The function returns a \code{data.frame} object with
    the installed local annotations.
}
\description{
    This function returns a data frame with information
    on locally installed, supported or custom, annotations.
    
}
\examples{
db <- file.path(system.file(package="metaseqR2"),
  "annotation.sqlite")
if (file.exists(db))
  ig <- getInstalledAnnotations(obj=db)
}
\author{
    Panagiotis Moulos
}

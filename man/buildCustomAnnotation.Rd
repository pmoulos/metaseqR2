\name{buildCustomAnnotation}
\alias{buildCustomAnnotation}
\title{Import custom annotation to the metaseqR2 annotation
    database from GTF file}
\usage{
    buildCustomAnnotation(gtfFile, metadata,
	db = file.path(system.file(package = "metaseqR2"),
	    "annotation.sqlite"), rewrite=TRUE)
}
\arguments{
    \item{gtfFile}{a GTF file containing the gene structure
    of the organism to be imported.}

    \item{metadata}{a list with additional information about
    the annotation to be imported. See Details.}
    
    \item{db}{a valid path (accessible at least by the
    current user) where the annotation database will be 
    set up. It defaults to 
    \code{system.file(package = "metaseqR2"),}
    \code{"annotation.sqlite")} that is, the installation
    path of metaseqR2 package. See also Details.}
    
    \item{rewrite}{if custom annotation found, rwrite? 
    (default \code{FALSE}). Set to \code{TRUE} if you wish to 
    update the annotation database for a particular custom
    annotation.}
}
\value{
    The function does not return anything. Only the SQLite 
    database is created or updated.
}
\description{
    This function imports a GTF file with some custom annotation
    to the metaseqR2 annotation database.
}
\details{
	Regarding the \code{metadata} argument, it is a list
	with specific format which instructs 
	\code{buildCustomAnnotation} on importing the custom
	annotation. Such a list may has the following members: 
	\itemize{
		\item \code{organism} a name of the organism which is
		imported (e.g. \code{"my_mm9"}). This is the only
		mandatory member.
		\item \code{source} a name of the source for this
		custom annotation (e.g. \code{"my_mouse_db"}). If
		not given or \code{NULL}, the word \code{"inhouse"}
		is used.
		\item \code{version} a string denoting the version.
		If not given or \code{NULL}, current date is used.
		\item \code{chromInfo} it can be one of the following:
		\itemize{
			\item a tab-delimited file with two columns, the
			first being the chromosome/sequence names and the
			second being the chromosome/sequence lengths.
			\item a BAM file to read the header from and
			obtain the required information
			\item a \code{\link{data.frame}} with one column
			with chromosome lengths and chromosome names as
			\code{rownames}.
		}
	}
	See the examples below for a \code{metadata} example.
	
	Regarding \code{db}, this controls the location of the
	installation database. If the default is used, then there is
	no need to provide the local database path to any function
	that uses the database (e.g. the main \code{metaseqr2}).
	Otherwise, the user will either have to provide this each
	time, or the annotation will have to be downloaded and used
	on-the-fly.
}
\examples{
\donttest{
    # FIXME: Example of custom annotation
    
    # Build on the default location
    buildCustomAnnotation(gtfFile, metadata,
		db = file.path(path.expand("~"), "my_ann.sqlite"))
}
}
\author{
    Panagiotis Moulos
}

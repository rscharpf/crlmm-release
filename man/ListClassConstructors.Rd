\name{ListClassConstructors}
\alias{BafLrrSetList}
\alias{OligoSetList}
\title{ Methods for Function \code{BafLrrSetList} in Package \pkg{crlmm} ~~}
\description{
  Constructors for \code{BafLrrSetList} and \code{OligoSetList} objects.
}

\usage{
BafLrrSetList(object, ...)
OligoSetList(object, ...)
}

\arguments{
  \item{object}{ A \code{CNSet} object.}
  \item{...}{Additional arguments \code{batch.name} and \code{chrom} can
  be used to specify specific batches or chromosomes in the \code{CNSet}
  object.
}
}

\details{

  Constructs a \code{BafLrrSetList} object or a \code{OligoSetList}
  object from an object of class \code{CNSet}.

}

\value{A \code{BafLrrSetList} or \code{OligoSetList}}

\examples{
data(cnSetExample)
oligoList <- OligoSetList(cnSetExample)
## only contains 1 chromosome, so list only has one element
dims(oligoList)
brList <- BafLrrSetList(cnSetExample)
dims(brList)
}

\seealso{\code{\linkS4class{BeadStudioSetList}}}


\keyword{manip}
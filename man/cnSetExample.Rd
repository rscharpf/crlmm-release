\name{cnSetExample}
\alias{cnSetExample}
\alias{cnSetExample2}
\docType{data}
\title{
	Object of class 'CNSet'
}
\description{

	The data for the first 16 polymorphic markers in the HapMap analysis.

}
\details{
  This object was created from the copynumber vignette in
  inst/scripts. A subset of markers was selected to keep the package
  size small.
}
\usage{

data(cnSetExample)
data(cnSetExample2)

}
\format{
  The data illustrates the \code{CNSet-class}, with
	\code{assayData} containing the quantile-normalized
	intensities for the A and B alleles, genotype calls and
	confidence scores.

}
\examples{
     data(cnSetExample)
     data(cnSetExample2)
}
\keyword{datasets}

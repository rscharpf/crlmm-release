\name{CNSet-methods}
\Rdversion{1.1}
\docType{methods}
\alias{coerce,CNSet,oligoSnpSet-method}
\alias{CA,CNSet-method}
\alias{CB,CNSet-method}
\alias{lines,CNSet-method}
\alias{totalCopynumber,CNSet-method}
\alias{rawCopynumber,CNSet-method}
\alias{nuA,CNSet-method}
\alias{nuB,CNSet-method}
\alias{Ns,CNSet-method}
\alias{corr,CNSet-method}
\alias{mads,CNSet-method}
\alias{medians,CNSet-method}
\alias{phiA,CNSet-method}
\alias{phiB,CNSet-method}
\alias{tau2,CNSet-method}
\alias{OligoSetList,CNSet-method}
\alias{BafLrrSetList,CNSet-method}

\title{crlmm methods for class "CNSet"}
\description{

	CNSet is a container defined in the oligoClasses package for
	 storing normalized intensities for genotyping platforms,
	 genotype calls, and parameters estimated for copy
	 number. Accessors for data that an object of this class
	 contains are largely defined in the package
	 oligoClasses. CNSet methods that involve more complex
	 calculations that are specific to the crlmm package, such as
	 computing allele-specific copy number, are included in crlmm
	 and described here.

}

\section{Methods}{
  \describe{
    \item{}{\code{as(from, "oligoSnpSet")}: Method for coercing object
      \code{from} (class \code{CNSet}) to an object of class
      \code{oligoSnpSet}.}
    \item{CA}{\code{signature(object="CNSet")}: calculates raw copy number
      for allele A}
    \item{CB}{\code{signature(object="CNSet")}: calculates raw copy
  number for allele B}
    \item{lines}{\code{signature(x="CNSet")}: plot ellipses (95th
  percentile) for prediction regions}
    \item{totalCopynumber}{\code{signature(object="CNSet")}: calculates
  total raw copy number}
    \item{rawCopynumber}{\code{signature(object="CNSet")}: same as totalCopynumber}
    \item{nuA}{\code{signature(object="CNSet")}: estimate of mean
  background (intensity-scale) for allele A}
    \item{nuB}{\code{signature(object="CNSet")}: estimate of mean
  background (intensity-scale) for allele A}
    \item{phiA}{\code{signature(object="CNSet")}: estimate of slope coefficient (intensity-scale) for allele A}
    \item{phiB}{\code{signature(object="CNSet")}: estimate of slope coefficient (intensity-scale) for allele B}
    \item{Ns}{\code{signature(object="CNSet")}: genotype frequencies }
    \item{corr}{\code{signature(object="CNSet")}: correlation of
  log-transformed normalized intensities within the genotype clusters}
    \item{mads}{\code{signature(x="CNSet")}: ...}
    \item{medians}{\code{signature(object="CNSet")}: ...}
    \item{tau2}{\code{signature(object="CNSet")}: ...}


  \item{}{\code{OligoSetList(object)}: constructs an object of class
    \code{OligoSetList} from \code{object} having class \code{CNSet}.
  }

  \item{}{\code{BafLrrSetList(object)}: constructs an object of class
    \code{BafLrrSetList} from \code{object} having class \code{CNSet}.
  }

  }
}

\seealso{
  \code{\link{CNSet-class}}, \code{\link{CA}}, \code{\link{CB}},
  \code{\link{totalCopynumber}}, \code{\link{rawCopynumber}}
}
\keyword{methods}

\name{AssayData-methods}
\Rdversion{1.1}
\docType{methods}
\alias{Ns,AssayData-method}
\alias{corr,AssayData-method}
\alias{mads,AssayData-method}
\alias{medians,AssayData-method}
\alias{tau2,AssayData-method}

\title{Methods for class "AssayData" in crlmm}

\description{

	The \code{batchStatistics} slot in a \code{CNSet} object is an
	instance of the \code{AssayData} slot. In general, the
	accessors for \code{AssayData} are called indirectly by the
	corresponding method for the \code{CNSet} class and not called
	directly by the user.

}

\section{Methods}{
  \describe{
    \item{Ns}{\code{signature(object="AssayData")}:  Accessor for genotype
  frequencies}
    \item{corr}{\code{signature(object="AssayData")}: Accessor for the correlation of the
  log-transformed normalized intensities within the diallelic genotype clusters}
    \item{mads}{\code{signature(x="AssayData")}: Accessor for the median
  absolute deviation of the  normalized intensities within the
  diallelic genotype clusters }
    \item{medians}{\code{signature(object="AssayData")}: Accessor for
  the posterior mean of the normalized intensity within the diallelic genotype
  clusters.}
    \item{tau2}{\code{signature(object="AssayData")}: Accessor for the
  median absolute deviation of the log-transformed intensities within
  the diallelic genotype clusters}
  }
}

\seealso{
  \code{\link{CNSet-class}}, \code{\link{Ns}}, \code{\link{tau2}}, \code{\link{corr}}, \code{\link{mads}}, \code{\link{medians}}
}
\keyword{manip}

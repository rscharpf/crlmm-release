\name{genotypeAffy}
\alias{genotypeAffy}
\title{Genotype Affymetrix CEL files}
\description{Assign diallelic genotypes at polymorphic markers}
\usage{
genotypeAffy(cnSet, SNRMin = 5, recallMin = 10, recallRegMin = 1000, gender = NULL, badSNP = 0.7, returnParams = TRUE, verbose = TRUE)
}

\arguments{
  \item{cnSet}{
    An object of class \code{CNSet}
  }
  \item{SNRMin}{
    See \code{\link{crlmm}}
  }
  \item{recallMin}{
    See \code{\link{crlmm}}
}
  \item{recallRegMin}{
    See \code{\link{crlmm}}
}
  \item{gender}{
    See \code{\link{crlmm}}
}
  \item{badSNP}{
    See \code{\link{crlmm}}
}
  \item{returnParams}{
    See \code{\link{crlmm}}
}
  \item{verbose}{
    Logical.
}
}
\details{
  Wrapper for crlmm genotyping.
}
\value{
  Returns logical.  SNP genotypes and confidence scores are written to
  \code{ff_matrix} objects.
}

\author{
R.Scharpf
}
\seealso{
  \code{\link{crlmm}}, \code{\link{calls}}, \code{\link{confs}}
}
\keyword{classif}
\name{cnrmaAffy}
\alias{cnrmaAffy}
\title{quantile normalize nonpolymorphic markers}
\description{Quantile normalize nonpolymorphic markers to hapmap reference distribution}
\usage{
cnrmaAffy(cnSet, seed = 1, verbose = TRUE)
}

\arguments{
  \item{cnSet}{
    Object of class \code{CNSet}
  }
  \item{seed}{
    Random number seed
  }
  \item{verbose}{
    Logical.
  }
}
\value{
  Returns logical.  Normalized intensities are written to the alleleA
  \code{ff_matrix} stored in the \code{CNSet} assayData.
}
\author{ R. Scharpf}

\seealso{\code{\link{snprmaAffy}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{robust}
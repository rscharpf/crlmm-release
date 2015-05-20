\name{snprmaAffy}
\alias{snprmaAffy}
\title{Quantile normalize intensities for SNPs}
\description{Quantile normalize intensities for SNPs to a HapMap target
  reference distribution}
\usage{
snprmaAffy(cnSet, mixtureSampleSize = 10^5, eps = 0.1, seed = 1, verbose = TRUE)
}

\arguments{
  \item{cnSet}{Object of class \code{CNSet}}
  \item{mixtureSampleSize}{Sample size to be use when fitting the
  mixture model. }
  \item{eps}{Stop criteria.}
  \item{seed}{Seed to be used when sampling. }
  \item{verbose}{Logical.}
}

\value{
  Returns nothing.  Normalized intensities are written to files.
}
\author{
R.Scharpf
}

\seealso{
  \code{\link{snprma}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{robust}

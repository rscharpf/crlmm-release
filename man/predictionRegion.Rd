\name{predictionRegion}
\alias{predictionRegion}
\alias{predictionRegion,CNSet,integer-method}
\title{Prediction regions for integer copy number}
\description{Bivariate normal prediction regions for integer copy
  number.  Copy numbers 0-4 allowed.}
\usage{
predictionRegion(object, copyNumber)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{CNSet} object.}
  \item{copyNumber}{Integer vector.  0-4 allowed.}
}
\details{

  We fit a linear regression for each allele to the diallic genotype
  cluster medians.  Denoting the background and slope by nu and phi,
  respectively, the mean for the bivariate normal prediction region is
  given by

  mu_A = nu_A + CA * phi_A

  and

  mu_B nu_B + CB * phi_B

  The variance and correlation of the normalized intensities is
  estimated from the diallelic genotype clusters AA, AB, and BB on the
  log-scale.  For copy number not equal to two, we assume that the
  variance is approximately the same for copy number not equal to 2.

}

\value{

  A list named by the genotype.  `NULL' refers to copy number zero, `A'
  is a hemizygous deletion, etc. Each element is a list of the means
  (mu) and covariance (cov) for each marker stored as an array.  For
  `mu', the dimensions of the array are marker x allele (A or B) x
  batch.  For `cov', the dimensions of the array are marker x 3
  (varA, cor, and varB) x batch.

}


\references{
  Scharpf et al., 2011, Biostatistics.
}
\author{
  R. Scharpf
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%  \code{\link{calculatePosteriorMean}},
  \code{\link{posteriorProbability}}, \code{\link{genotypes}}
}
\examples{
data(cnSetExample)
pr <- predictionRegion(cnSetExample, copyNumber=0:4)
names(pr)
## bivariate normal prediction region for NULL genotype (homozygous deletion)
str(pr[["NULL"]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{list}
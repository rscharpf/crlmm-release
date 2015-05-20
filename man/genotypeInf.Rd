\name{genotypeInf}
\alias{genotypeInf}

\title{

	Genotyping of Illumina Infinium II arrays.

}
\description{

	Genotyping of Illumina Infinium II arrays. This function
	provides CRLMM/KRLMM genotypes and confidence scores for the the
	polymorphic markers and is a required step prior to copy
	number estimation.

}
\usage{
genotypeInf(cnSet, mixtureParams, probs = rep(1/3, 3), SNRMin = 5,
recallMin = 10, recallRegMin = 1000, verbose = TRUE, returnParams =
TRUE, badSNP = 0.7, gender = NULL, DF = 6, cdfName, call.method="crlmm", 
trueCalls = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cnSet}{An object of class \code{CNSet}}
  \item{mixtureParams}{
    \code{data.frame} containing mixture model parameters needed for
    genotyping.  The mixture model parameters are estimated from the
    \code{preprocessInf} function.
  }
  \item{probs}{'numeric' vector with priors for AA, AB and BB.}
  \item{SNRMin}{'numeric' scalar defining the minimum SNR used to filter
    out samples.}
  \item{recallMin}{Minimum number of samples for recalibration. }
  \item{recallRegMin}{Minimum number of SNP's for regression.}
  \item{verbose}{  'logical.'  Whether to print descriptive messages during processing.}
  \item{returnParams}{'logical'. Return recalibrated parameters from crlmm.}
  \item{badSNP}{'numeric'. Threshold to flag as bad SNP (affects
  batchQC)}
  \item{gender}{  integer vector (  male = 1, female =2 ) or missing,
  with same length as filenames.  If missing, the gender is
  predicted.}
  \item{DF}{'integer' with number of degrees of freedom to use with
    t-distribution.}
  \item{cdfName}{\code{character} string indicating which annotation
    package to load.}
  \item{call.method}{character string specifying the genotype calling algorithm to use ('crlmm' or 'krlmm').}
  \item{trueCalls}{ matrix specifying known Genotype calls for a subset of samples and features(1 - AA, 2 - AB, 3 - BB).}
}

\details{

	The genotype calls and confidence scores are written to
	file using \code{ff} protocols for I/O.  For the most part,
	the calls and confidence scores can be accessed as though the
	data is in memory through the methods \code{snpCall} and
	\code{snpCallProbability}, respectively.

	The genotype calls are stored using an integer representation: 1
	- AA, 2 - AB, 3 - BB.  Similarly, the call probabilities are
	stored using an integer representation to reduce file size using
	the transformation 'round(-1000*log2(1-p))', where p is the
	probability.  The function \code{i2P} can be used to convert the
	integers back to the scale of probabilities.

	An optional \code{trueCalls} argument can be provided to KRLMM method
	which contains known genotype calls(can contain some NAs) for some 
	samples and SNPs. This will used to compute KRLMM parameters by calling
	\code{vglm} function from \code{VGAM} package.

	The KRLMM method makes use of functions provided in \code{parallel} 
	package to speed up the process. It by default initialises up to 
	8 clusters. This is configurable by setting up an option named 
	"krlmm.cores", e.g. options("krlmm.cores" = 16). 
      }

\value{

Logical.  If the genotyping is completed, the value 'TRUE' is
returned.  Note that \code{assayData} elements 'call' and
'callProbability' are updated on disk. Therefore, the genotypes and
confidence scores can be retrieved using accessors for the
\code{CNSet}  class.

}
\author{
R. Scharpf
}
\seealso{
	\code{\link{crlmm}}, \code{\link{snpCall}},
	\code{\link{snpCallProbability}}, \code{\link{annotationPackages}}
}
\examples{
	## See the 'illumina_copynumber' vignette in inst/scripts of
	## the source package
}
\keyword{classif}
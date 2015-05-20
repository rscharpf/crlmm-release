\name{crlmmCopynumber}
\alias{crlmmCopynumber}
\alias{crlmmCopynumber2}
\alias{crlmmCopynumberLD}
\title{Locus- and allele-specific estimation of copy number}
\description{
  Locus- and allele-specific estimation of copy number.
}
\usage{
crlmmCopynumber(object, MIN.SAMPLES=10, SNRMin = 5, MIN.OBS = 1,
	        DF.PRIOR = 50, bias.adj = FALSE,
                prior.prob = rep(1/4, 4), seed = 1, verbose = TRUE,
                GT.CONF.THR = 0.80, MIN.NU = 2^3, MIN.PHI = 2^3,
                THR.NU.PHI = TRUE, type=c("SNP", "NP", "X.SNP", "X.NP"),
                fit.linearModel=TRUE)
}

\arguments{

  \item{object}{object of class \code{CNSet}.}

 \item{MIN.SAMPLES}{ 'Integer'.  The minimum number of samples in a
  batch.  Bathes with fewer than MIN.SAMPLES are skipped.  Therefore,
  samples in batches with fewer than MIN.SAMPLES have NA's for the
  allele-specific copy number and NA's for the linear model
  parameters.

}

  \item{SNRMin}{ Samples with low signal to noise ratios are  excluded.  }

  \item{MIN.OBS}{

  For a SNP with with fewer than \code{MIN.OBS} of a genotype in a given
  batch, the within-genotype median is imputed.  The imputation is based
  on a regression using SNPs for which all three biallelic genotypes are
  observed.  For example, assume at at a given SNP genotypes AA and AB
  were observed and BB is an unobserved genotype.  For SNPs in which all
  3 genotypes were observed, we fit the model E(mean_BB) = beta0 +
  beta1*mean_AA + beta2*mean_AB, obtaining estimates; of beta0, beta1,
  and beta2.  The imputed mean at the SNP with unobserved BB is then
  beta0hat + beta1hat * mean_AA of beta2hat * mean_AB.

}

  \item{DF.PRIOR}{

  The 2 x 2 covariance matrix of the background and signal variances
  is estimated from the data at each locus.  This matrix is then
  smoothed towards a common matrix estimated from all of the loci.
  DF.PRIOR controls the amount of smoothing towards the common matrix,
  with higher values corresponding to greater smoothing.  Currently,
  DF.PRIOR is not estimated from the data.  Future versions may
  estimate DF.PRIOR empirically.

}

\item{bias.adj}{

  \code{bias.adj} is currently ignored (as well as the prior.prob
  argument).  We plan to add this feature back to the crlmm package in
  the near future. This feature, when \code{TRUE}, updated initial
  estimates from the linear model after excluding samples with a low
  posterior probability of normal copy number.  Excluding samples that
  have a low posterior probability can be helpful at loci in which a
  substantial fraction of the samples have a copy number alteration.
  For additional information, see Scharpf et al., 2010.

}
  \item{prior.prob}{

    This argument is currently ignored.  A numerical vector providing
  prior probabilities for copy number states corresponding to homozygous
  deletion, hemizygous deletion, normal copy number, and amplification,
  respectively.

}
  \item{seed}{ Seed for random number generation.}

  \item{verbose}{ Logical. }

  \item{GT.CONF.THR}{

    Confidence threshold for genotype calls (0, 1).  Calls with
    confidence scores below this theshold are not used to estimate the
    within-genotype medians. See Carvalho et al., 2007 for information
    regarding confidence scores of biallelic genotypes.

}


  \item{MIN.NU}{ numeric. Minimum value for background
    intensity. Ignored if \code{THR.NU.PHI} is \code{FALSE}. }

  \item{MIN.PHI}{numeric. Minimum value for slope. Ignored if
    \code{THR.NU.PHI} is \code{FALSE}.}

  \item{THR.NU.PHI}{  If \code{THR.NU.PHI} is \code{FALSE},
    \code{MIN.NU} and \code{MIN.PHI} are ignored. When TRUE, background
    (nu) and slope (phi) coefficients below MIN.NU and MIN.PHI are set
    to MIN.NU and MIN.PHI, respectively.}

  \item{type}{ Character string vector that must be one or more of
    "SNP", "NP", "X.SNP", or "X.NP". Type refers to a set of
    markers. See details below}

  \item{fit.linearModel}{Logical.  If TRUE, a linear model is fit to
  estimate the parameters for computing the absolute copy number.  If
  FALSE, we compute the batch-specific, within-genotype median and MAD
  at polymorphic loci and the median and MAD at nonpolymorphic loci.}

 }

 \references{

  Carvalho B, Bengtsson H, Speed TP, Irizarry RA. Exploration,
  normalization, and genotype calls of high-density oligonucleotide SNP
  array data. Biostatistics. 2007 Apr;8(2):485-99. Epub 2006 Dec
  22. PMID: 17189563.

  Carvalho BS, Louis TA, Irizarry RA.
  Quantifying uncertainty in genotype calls.
  Bioinformatics. 2010 Jan 15;26(2):242-9.

  Scharpf RB, Ruczinski I, Carvalho B, Doan B, Chakravarti A, and
  Irizarry RA, Biostatistics.  Biostatistics, Epub July 2010.

}

\details{

    We suggest a minimum of 10 samples per batch for using
    \code{crlmmCopynumber}.  50 or more samples per batch is preferred
    and will improve the estimates.

    The functions \code{crlmmCopynumberLD} and
    \code{crlmmCopynumber2} have been deprecated.

	The argument \code{type} can be used to specify a subset of
	markers for which the copy number estimation algorithm is run.
	One or more of the following possible entries are valid: 'SNP',
	'NP', 'X.SNP', and 'X.NP'.

	'SNP' referers to autosomal SNPs.

	'NP' refers to autosomal nonpolymorphic markers.

	'X.SNP' refers to SNPs on chromosome X.

	'X.NP' refers to autosomes on chromosome X.

	However, users must run 'SNP' prior to running 'NP' and 'X.NP',
	or specify \code{type = c('SNP', 'X.NP')}.

      }

      \value{

	The value returned by the \code{crlmmCopynumber} function
	depends on whether the data is stored in RAM or whether the data
	is stored on disk using the R package \code{ff} for reading /
	writing.  If uncertain, the first line of the \code{show} method
	defined for \code{CNSet} objects prints whether the
	\code{assayData} elements are derived from the \code{ff} package
	in the first line.  Specifically,

	- if the elements of the \code{batchStaticts} slot in the
	\code{CNSet} object have the class "ff_matrix" or "ffdf", then
	the \code{crlmmCopynumber} function updates the data stored on
	disk and returns the value \code{TRUE}.

	- if the elements of the \code{batchStatistics} slot in the
	\code{CNSet} object have the class 'matrix', then the
	\code{crlmmCopynumber} function returns an object of class
	\code{CNSet} with the elements of \code{batchStatistics}
	updated.

      }

\author{R. Scharpf}
\keyword{manip}

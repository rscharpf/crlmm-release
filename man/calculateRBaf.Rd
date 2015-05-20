\name{calculateRBaf}
\alias{calculateRBaf}
\alias{calculateRBaf,CNSet-method}
\title{Calculate log R ratios and B allele frequencies.}
\description{Calculate log R ratios and B allele frequencies from
  a \code{CNSet} object}
\usage{
calculateRBaf(object, batch.name, chrom)
}

\arguments{

  \item{object}{A \code{CNSet} object.}

  \item{batch.name}{A character string indicating the batch.  If
  missing, log R ratios and B allele frequencies are calculated for all
  batches in the \code{object}.}

  \item{chrom}{Integer indicating which chromosome to process. If
  missing, B allele frequencies and log R ratios are calculated for all
  autosomal chromosomes and chromosome X that are included in
  \code{object}.}

}

\details{

  \code{batch.name} must be a value in \code{batch(object)}. Currently,
  one must specify a single \code{batch.name}. If a character vector for
  \code{batch.name} is supplied, only the first is evaluated.

  TODO: A description of how these values are calculated.

}

\value{
  A named list.

  \code{baf}: Each element in the baf list is a matrix of B allele
  frequencies (one matrix for each chromosome).

  \code{lrr}: Each element in the lrr list is a matrix of log R ratios
  (one matrix for each chromosome).

  The log R ratios were scaled by a factor of 100 and stored as an
  integer.  B allele frequencies were scaled by a factor of 1000 and
  stored as an integer.

}

\references{

  Peiffer et al., High-resolution genomic profiling of chromosomal aberrations using
  Infinium whole-genome genotyping (2006), Genome Research

}

\author{Lynn Mireless}

\examples{
	data(cnSetExample)
	baf.lrr <- suppressWarnings(calculateRBaf(cnSetExample, "SHELF"))
	hist(baf.lrr[["baf"]][[1]]/1000, breaks=100)
	hist(baf.lrr[["lrr"]][[1]]/100, breaks=100)
	\dontrun{
		library(ff)
		baf.lrr <- suppressWarnings(calculateRBaf(cnSetExample, "SHELF"))
		class(baf.lrr[["baf"]][[1]]) ## ff_matrix
		class(baf.lrr[["lrr"]][[1]]) ## ff_matrix
	}
}
\keyword{list}
\keyword{methods}

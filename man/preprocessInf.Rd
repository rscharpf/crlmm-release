\name{preprocessInf}
\alias{preprocessInf}
\title{

	Preprocessing of Illumina Infinium II arrays.

}
\description{

	This function normalizes the intensities for the 'A' and 'B'
	alleles for a \code{CNSet} object and estimates mixture
	parameters used for subsequent genotyping.  See details for
	how the normalized intensities are written to file.  This step
	is required for subsequent genotyping and copy number
	estimation.

}
\usage{
preprocessInf(cnSet, sampleSheet=NULL, arrayNames = NULL, ids = NULL,
path = ".", arrayInfoColNames = list(barcode = "SentrixBarcode_A",
position = "SentrixPosition_A"), highDensity = TRUE, sep = "_", fileExt
= list(green = "Grn.idat", red = "Red.idat"), XY, saveDate = TRUE, stripNorm
= TRUE, useTarget = TRUE, mixtureSampleSize = 10^5, fitMixture = TRUE, 
quantile.method="between", eps = 0.1, verbose = TRUE, seed = 1, cdfName)
}

\arguments{
  \item{cnSet}{
  object of class \code{CNSet}
}

  \item{sampleSheet}{\code{data.frame} containing Illumina sample sheet
    information (for required columns, refer to BeadStudio Genotyping
    guide - Appendix A).}

  \item{arrayNames}{character vector containing names of arrays to be
    read in.  If \code{NULL}, all arrays that can be found in the
    specified working directory will be read in.}

  \item{ids}{vector containing ids of probes to be read in.  If
    \code{NULL} all probes found on the first array are read in.}

  \item{path}{character string specifying the location of files to be
    read by the function}

  \item{arrayInfoColNames}{(used when \code{sampleSheet} is specified)
    list containing elements 'barcode' which indicates column names in
    the \code{sampleSheet} which contains the arrayNumber/barcode number
    and 'position' which indicates the strip number.  In older style
    sample sheets, this information is combined (usually in a column
    named 'SentrixPosition') and this should be specified as
    \code{list(barcode=NULL, position="SentrixPosition")}}

  \item{highDensity}{logical (used when \code{sampleSheet} is
    specified). If \code{TRUE}, array extensions '\_A', '\_B' in
    sampleSheet are replaced with 'R01C01', 'R01C02' etc.}
  \item{sep}{character string specifying separator used in .idat file
    names.}
  \item{fileExt}{list containing elements 'Green' and 'Red' which
    specify the .idat file extension for the Cy3 and Cy5 channels.}
  \item{XY}{an \code{NChannelSet} object containing X and Y intensities.}
  \item{saveDate}{'logical'.  Should the dates from each .idat be saved
    with sample information?}
  \item{stripNorm}{'logical'.  Should the data be strip-level normalized?}
  \item{useTarget}{'logical' (only used when \code{stripNorm=TRUE}).
    Should the reference HapMap intensities be used in strip-level normalization?}
  \item{mixtureSampleSize}{ Sample size to be use when fitting the mixture model.}
  \item{fitMixture}{ 'logical.' Whether to fit per-array mixture
    model.}
  \item{quantile.method}{character string specifying the quantile normalization method to use ('within' or 'between' channels).}
  \item{eps}{   Stop criteria.}
  \item{verbose}{  'logical.'  Whether to print descriptive messages during processing.}
  \item{seed}{ Seed to be used when sampling. Useful for
    reproducibility}
  \item{cdfName}{ \code{character} string indicating which annotation
    package to load.}
}
\details{

	The normalized intensities are written to disk using package
	\code{ff} protocols for writing/reading to disk. Note that the
	object \code{CNSet} containing the \code{ff} objects in the
	\code{assayData} slot will be updated after applying this
	function.

}
\value{

	A \code{ff_matrix} object containing parameters for fitting the
	mixture model.  Note that while the \code{CNSet} object is not
	returned by this function, the object will be updated as the
	normalized intensities are written to disk.  In particular,
	after applying this function the normalized intensities in the
	\code{alleleA} and \code{alleleB} elements of \code{assayData}
	are now available.

}
\author{
R. Scharpf
}
\seealso{
	\code{\link{CNSet-class}}, \code{\link{A}}, \code{\link{B}},
	\code{\link{constructInf}}, \code{\link{genotypeInf}}, \code{\link{annotationPackages}}
}
\examples{
	## See the 'illumina_copynumber' vignette in inst/scripts of
	## the source package
}
\keyword{manip}
\name{constructInf}
\alias{constructInf}
\title{
	Instantiate an object of class CNSet for the Infinium platforms.
}
\description{
	Instantiates an object of class CNSet for the Infinium
	platforms.  Elements of \code{assayData} and
	\code{batchStatistics} will be \code{ff} objects. See details.
}
\usage{
constructInf(sampleSheet = NULL, arrayNames = NULL, path = ".",
       arrayInfoColNames = list(barcode="SentrixBarcode_A",position="SentrixPosition_A"), highDensity = FALSE, sep = "_",
       fileExt = list(green = "Grn.idat", red = "Red.idat"), XY, cdfName, verbose = FALSE, batch=NULL, saveDate = TRUE)
}
\arguments{
  \item{sampleSheet}{\code{data.frame} containing Illumina sample sheet
    information (for required columns, refer to BeadStudio Genotyping
    guide - Appendix A).
}
  \item{arrayNames}{character vector containing names of arrays to be
    read in.  If \code{NULL}, all arrays that can be found in the
    specified working directory will be read in.}
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
  \item{XY}{an \code{NChannelSet} containing X and Y intensities.}
  \item{cdfName}{ annotation package  (see also \code{validCdfNames})}
  \item{verbose}{  'logical.'  Whether to print descriptive messages
  during processing.}
 \item{batch}{ batch variable. See details.}
  \item{saveDate}{'logical'.  Should the dates from each .idat be saved
    with sample information?}
}
\details{

	This function initializes a container for storing the normalized
	intensities for the A and B alleles at polymorphic loci and the
	normalized intensities for the 'A' allele at nonpolymorphic
	loci. CRLMM genotype calls and confidence scores are also stored
	in the assayData. This function does not do any preprocessing or
	genotyping -- it only creates an object of the appropriate size.
	The initialized values will all be 'NA'.

	The ff package provides infrastructure for accessing and writing
	data to disk instead of keeping data in memory.  Each element of
	the \code{assayData} and \code{batchStatistics} slot are ff
	objects.  ff objects in the R workspace contain pointers to
	several files with the '.ff' extension on disk.  The location of
	where the data is stored on disk can be specified by use of the
	\code{ldPath} function.  Users should not move or rename this
	directory.  If only output files are stored in \code{ldPath},
	one can either remove the entire directory prior to rerunning
	the analysis or all of the '.ff' files.  Otherwise, one would
	accumulate a large number of '.ff' files on disk that are no
	longer in use.

	We have adopted the \code{ff} package in order to reduce crlmm's
	memory footprint. The memory usage can be fine-tuned by the
	utilities \code{ocSamples} and \code{ocProbesets} provided in
	the \code{oligoClasses} package. In most instances, the
	user-level interface will be no different than accessing data
	from ordinary matrices in R.  However, the differences in the
	underlying representation can become more noticeable for very
	large datasets in which the I/O for accessing data from the disk
	can be substantial.

}
\value{

	A \code{CNSet} object

}
\author{
R. Scharpf
}
\seealso{
	\code{\link{ldPath}}, \code{\link{ocSamples}}, \code{\link{ocProbesets}}, \code{\link{CNSet-class}}, \code{\link{preprocessInf}}, \code{\link{genotypeInf}}
}
\examples{
## See the Illumina vignettes in inst/scripts of the
## source package for an example
}
\keyword{manip}

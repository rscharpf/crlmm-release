\name{crlmmIlluminaV2}
\alias{crlmmIlluminaV2}
\title{Read and Genotype Illumina Infinium II BeadChip data with CRLMM}
\description{
  Implementation of the CRLMM algorithm for
  data from Illumina's Infinium II BeadChips.
}
\usage{

crlmmIlluminaV2(sampleSheet=NULL, arrayNames=NULL, ids=NULL, path=".",
      arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
      highDensity=FALSE, sep="_", fileExt=list(green="Grn.idat", red="Red.idat"),
      saveDate=FALSE, stripNorm=TRUE, useTarget=TRUE,
      row.names=TRUE, col.names=TRUE, probs=c(1/3, 1/3, 1/3),
      DF=6, SNRMin=5, gender=NULL, seed=1, mixtureSampleSize=10^5,
      eps=0.1, verbose=TRUE, cdfName, sns, recallMin=10,
      recallRegMin=1000, returnParams=FALSE, badSNP=.7)
}

\arguments{
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
  \item{saveDate}{'logical'.  Should the dates from each .idat be saved
    with sample information?}
  \item{stripNorm}{'logical'.  Should the data be strip-level normalized?}
  \item{useTarget}{'logical' (only used when \code{stripNorm=TRUE}).
    Should the reference HapMap intensities be used in strip-level normalization?}
  \item{row.names}{'logical'. Use rownames - SNP names?}
  \item{col.names}{'logical'. Use colnames - Sample names?}
  \item{probs}{'numeric' vector with priors for AA, AB and BB.}
  \item{DF}{'integer' with number of degrees of freedom to use with t-distribution.}
  \item{SNRMin}{'numeric' scalar defining the minimum SNR used to filter
  out samples.}
  \item{gender}{'integer' vector, with same length as 'filenames',
    defining sex. (1 - male; 2 - female)}
  \item{seed}{'integer' scalar for random number generator (used to
    sample \code{mixtureSampleSize} SNPs for mixture model.}
  \item{mixtureSampleSize}{'integer'. The number of SNP's to be used
    when fitting the mixture model.}
  \item{eps}{Minimum change for mixture model.}
  \item{verbose}{'logical'.}
  \item{cdfName}{'character' defining the chip annotation (manifest) to use
    ('human370v1c', human550v3b', 'human650v3a', 'human1mv1c',
    'human370quadv3c', 'human610quadv1b', 'human660quadv1a',
    'human1mduov3b', 'humanomni1quadv1b', 'humanomniexpress12v1b', 'humancytosnp12v2p1h')}
  \item{sns}{'character' vector with sample names to be used.}
  \item{recallMin}{'integer'. Minimum number of samples for recalibration.}
  \item{recallRegMin}{'integer'. Minimum number of SNP's for regression.}
  \item{returnParams}{'logical'. Return recalibrated parameters.}
  \item{badSNP}{'numeric'. Threshold to flag as bad SNP (affects batchQC)}
}
\value{
  A \code{SnpSet} object which contains
  \item{calls}{Genotype calls (1 - AA, 2 - AB, 3 - BB)}
  \item{callProbability}{confidence scores 'round(-1000*log2(1-p))'}
  in the \code{assayData} slot and
  \item{SNPQC}{SNP Quality Scores}
  \item{batchQC}{Batch Quality Scores}
  along with center and scale parameters when \code{returnParams=TRUE}
  in the \code{featureData} slot.
}

\details{
  This function combines the reading of data from idat files using
  \code{readIdatFiles} and genotyping to reduce memory usage.
}

\references{
  Ritchie ME, Carvalho BS, Hetrick KN, Tavar\'{e} S, Irizarry RA.
  R/Bioconductor software for Illumina's Infinium whole-genome
  genotyping BeadChips. Bioinformatics. 2009 Oct 1;25(19):2621-3.

  Carvalho B, Bengtsson H, Speed TP, Irizarry RA. Exploration,
  normalization, and genotype calls of high-density oligonucleotide SNP
  array data. Biostatistics. 2007 Apr;8(2):485-99. Epub 2006 Dec
  22. PMID: 17189563.

  Carvalho BS, Louis TA, Irizarry RA.
  Quantifying uncertainty in genotype calls.
  Bioinformatics. 2010 Jan 15;26(2):242-9.
}

\author{Matt Ritchie}

\examples{
## crlmmOut = crlmmIlluminaV2(samples,path=path,arrayInfoColNames=list(barcode="Chip",position="Section"),
##                             saveDate=TRUE,cdfName="human370v1c",returnParams=TRUE)

}
\seealso{\code{\link{crlmmIllumina}}}
\keyword{classif}

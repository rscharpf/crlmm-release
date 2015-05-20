\name{crlmm}
\alias{crlmm}
\alias{crlmm2}

\title{Genotype oligonucleotide arrays with CRLMM}
\description{
  This is a faster and more efficient implementation of the CRLMM
  algorithm, especially designed for Affymetrix SNP 5 and 6 arrays (to
  be soon extended to other platforms).
}
\usage{
crlmm(filenames, row.names=TRUE, col.names=TRUE,
      probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5,
      gender=NULL, save.it=FALSE, load.it=FALSE,
      intensityFile, mixtureSampleSize=10^5,
      eps=0.1, verbose=TRUE, cdfName, sns, recallMin=10,
      recallRegMin=1000, returnParams=FALSE, badSNP=0.7)
crlmm2(filenames, row.names=TRUE, col.names=TRUE,
      probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5,
      gender=NULL, save.it=FALSE, load.it=FALSE,
      intensityFile, mixtureSampleSize=10^5,
      eps=0.1, verbose=TRUE, cdfName, sns, recallMin=10,
      recallRegMin=1000, returnParams=FALSE, badSNP=0.7)
}

\arguments{
  \item{filenames}{'character' vector with CEL files to be genotyped.}
  \item{row.names}{'logical'. Use rownames - SNP names?}
  \item{col.names}{'logical'. Use colnames - Sample names?}
  \item{probs}{'numeric' vector with priors for AA, AB and BB.}
  \item{DF}{'integer' with number of degrees of freedom to use with t-distribution.}
  \item{SNRMin}{'numeric' scalar defining the minimum SNR used to filter
  out samples.}
  \item{gender}{'integer' vector, with same length as 'filenames',
    defining sex. (1 - male; 2 - female)}
  \item{save.it}{'logical'. Save preprocessed data?}
  \item{load.it}{'logical'. Load preprocessed data to speed up analysis?}
  \item{intensityFile}{'character' with filename to be saved/loaded -
    preprocessed data.}
  \item{mixtureSampleSize}{Number of SNP's to be used with the mixture model.}
  \item{eps}{Minimum change for mixture model.}
  \item{verbose}{'logical'.}
  \item{cdfName}{'character' defining the CDF name to use
    ('GenomeWideSnp5', 'GenomeWideSnp6')}
  \item{sns}{'character' vector with sample names to be used.}
  \item{recallMin}{Minimum number of samples for recalibration.}
  \item{recallRegMin}{Minimum number of SNP's for regression.}
  \item{returnParams}{'logical'. Return recalibrated parameters.}
  \item{badSNP}{'numeric'. Threshold to flag as bad SNP (affects batchQC)}
}
\value{
  A \code{SnpSet} object.
  \item{calls}{Genotype calls (1 - AA, 2 - AB, 3 - BB)}
  \item{confs}{Confidence scores 'round(-1000*log2(1-p))'}
  \item{SNPQC}{SNP Quality Scores}
  \item{batchQC}{Batch Quality Score}
  \item{params}{Recalibrated parameters}
}
\details{
  'crlmm2' allows one to genotype very large datasets (via ff package) and also permits
  the use of clusters or multiple cores (via snow package) to speed up
  genotyping.

  As noted above, the call probabilities are stored using an integer
  representation to reduce file size using the transformation
  'round(-1000*log2(1-p))', where p is the probability.  The function
  \code{i2P} can be used to convert the integers back to the scale of
  probabilities.
}


\seealso{ \code{\link{i2p}}, \code{\link{snpCall}}, \code{\link{snpCallProbability}}}

\references{ Carvalho B, Bengtsson H, Speed TP,
  Irizarry RA. Exploration, normalization, and genotype calls of
  high-density oligonucleotide SNP array data. Biostatistics. 2007
  Apr;8(2):485-99. Epub 2006 Dec 22. PMID: 17189563.

  Carvalho BS, Louis TA, Irizarry RA.
  Quantifying uncertainty in genotype calls.
  Bioinformatics. 2010 Jan 15;26(2):242-9.
}
\examples{
## this can be slow
library(oligoClasses)
if (require(genomewidesnp6Crlmm) & require(hapmapsnp6)){
  path <- system.file("celFiles", package="hapmapsnp6")

  ## the filenames with full path...
  ## very useful when genotyping samples not in the working directory
  cels <- list.celfiles(path, full.names=TRUE)
  (crlmmOutput <- crlmm(cels))
  ## If gender is known, one should check that the assigned gender is
  ## correct, or pass the integer coding of gender as an argument to the
  ## crlmm function as done below
}

\dontrun{
## HPC Example
library(ff)
library(snow)
library(crlmm)
## genotype 50K SNPs at a time
ocProbesets(50000)
## setup cluster - 8 cores on the machine
library(doSNOW)
cl <- makeCluster(8, "SOCK")
registerDoSNOW(cl)
##setCluster(8, "SOCK")

path <- system.file("celFiles", package="hapmapsnp6")
cels <- list.celfiles(path, full.names=TRUE)
crlmmOutput <- crlmm2(cels)
}

}
\keyword{classif}


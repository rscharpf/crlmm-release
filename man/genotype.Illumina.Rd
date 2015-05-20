\name{genotype.Illumina}
\alias{genotype.Illumina}

\title{
	Preprocessing and genotyping of Illumina Infinium II arrays.
}
\description{
	Preprocessing and genotyping of Illumina Infinium II arrays.
}
\usage{
genotype.Illumina(sampleSheet=NULL, arrayNames=NULL, ids=NULL, path=".",
      arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
      highDensity=FALSE, sep="_", fileExt=list(green="Grn.idat", red="Red.idat"), XY=NULL,
      call.method="crlmm", trueCalls=NULL, cdfName, copynumber=TRUE, batch=NULL, saveDate=FALSE, stripNorm=TRUE, 
      useTarget=TRUE, quantile.method="between", mixtureSampleSize=10^5, fitMixture=TRUE,                               
      eps =0.1, verbose = TRUE, seed = 1, sns, probs = rep(1/3, 3), DF = 6, SNRMin = 5, 
      recallMin = 10, recallRegMin = 1000, gender = NULL, returnParams = TRUE, badSNP = 0.7)
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
  \item{XY}{\code{NChannelSet} containing X and Y intensities.}
  \item{call.method}{character string specifying the genotype calling algorithm to use ('crlmm' or 'krlmm').}
  \item{trueCalls}{matrix specifying known Genotype calls(can contain some NAs) for a subset of samples and features (1 - AA, 2 - AB, 3 - BB).}
  \item{cdfName}{ annotation package  (see also \code{validCdfNames})}
  \item{copynumber}{ 'logical.' Whether to store copy number intensities with SNP output.}
  \item{batch}{ character vector indicating the batch variable. Must be
	the same length as the number of samples. See details.}
  \item{saveDate}{'logical'.  Should the dates from each .idat be saved
    with sample information?}
  \item{stripNorm}{'logical'.  Should the data be strip-level normalized?}
  \item{useTarget}{'logical' (only used when \code{stripNorm=TRUE}).
    Should the reference HapMap intensities be used in strip-level normalization?}
  \item{quantile.method}{character string specifying the quantile normalization method to use ('within' or 'between' channels).}
  \item{mixtureSampleSize}{ Sample size to be use when fitting the mixture model.}
  \item{fitMixture}{ 'logical.' Whether to fit per-array mixture model.}
  \item{eps}{   Stop criteria.}
  \item{verbose}{  'logical.'  Whether to print descriptive messages during processing.}
  \item{seed}{ Seed to be used when sampling. Useful for reproducibility}
  \item{sns}{The sample identifiers.  If missing, the default sample names are \code{basename(filenames)}}
  \item{probs}{'numeric' vector with priors for AA, AB and BB.}
  \item{DF}{'integer' with number of degrees of freedom to use with t-distribution.}
  \item{SNRMin}{'numeric' scalar defining the minimum SNR used to filter
  out samples.}
  \item{recallMin}{Minimum number of samples for recalibration. }
  \item{recallRegMin}{Minimum number of SNP's for regression.}
  \item{gender}{  integer vector (  male = 1, female = 2 ) or missing,
  with same length as filenames.  If missing, the gender is predicted.}
  \item{returnParams}{'logical'. Return recalibrated parameters from crlmm.}
  \item{badSNP}{'numeric'. Threshold to flag as bad SNP (affects batchQC)}
}

\details{

	For large datasets it is important to utilize the large data
	support by installing and loading the ff package before calling
	the \code{genotype} function. In previous versions of the
	\code{crlmm} package, we used different functions for
	genotyping depending on whether the ff package is loaded, namely
	\code{genotype} and \code{genotype2}.  The \code{genotype}
	function now handles both instances.

	\code{genotype.Illumina} is a wrapper of the \code{crlmm}
	function for genotyping.  Differences include (1) that the copy
	number probes (if present) are also quantile-normalized and (2)
	the class of object returned by this function, \code{CNSet}, is
	needed for subsequent copy number estimation.  Note that the
	batch variable (a character string) that must be passed to this
	function has no effect on the normalization or genotyping steps.
	Rather, \code{batch} is required in order to initialize a
	\code{CNSet} container with the appropriate dimensions.

        The new 'krlmm' option is available for certain chip types. Optional 
	argument \code{trueCalls} matrix contains known Genotype calls 
	(1 - AA, 2 - AB, 3 - BB) for a subset of samples and features. This 
	will used to compute KRLMM coefficients by calling \code{vglm} function 
	from \code{VGAM} package.

	The 'krlmm' method makes use of functions provided in \code{parallel} 
	package to speed up the process. It by default initialises up to 
	8 clusters. This is configurable by setting up an option named 
	"krlmm.cores", e.g. options("krlmm.cores" = 16). 

      }

\value{	A \code{SnpSuperSet} instance.}
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

\author{Matt Ritchie, Cynthia Liu, Zhiyin Dai}

  \note{For large datasets, load the 'ff' package prior to genotyping
-- this will greatly reduce the RAM required for big jobs.  See
\code{ldPath} and \code{ocSamples}.  The function
\code{genotype.Illumina} supports parallelization, as the (not run)
example below indicates.}

\seealso{
	\code{\link{crlmmIlluminaV2}},
	\code{\link[oligoClasses]{ocSamples}},
	\code{\link[oligoClasses]{ldOpts}}
}
\examples{
\dontrun{
	# example for 'crlmm' option
	library(ff)
	library(crlmm)
	## to enable paralellization, set to TRUE
	if(FALSE){
		library(snow)
		library(doSNOW)
		## with 10 workers
		cl <- makeCluster(10, type="SOCK")
		registerDoSNOW(cl)
	}
	## path to idat files
	datadir <- "/thumper/ctsa/snpmicroarray/illumina/IDATS/370k"
	## read in your samplesheet
	samplesheet = read.csv(file.path(datadir, "HumanHap370Duo_Sample_Map.csv"), header=TRUE, as.is=TRUE)
	samplesheet <- samplesheet[-c(28:46,61:75,78:79), ]
	arrayNames <- file.path(datadir, unique(samplesheet[, "SentrixPosition"]))
	arrayInfo <- list(barcode=NULL, position="SentrixPosition")
	cnSet <- genotype.Illumina(sampleSheet=samplesheet,
				   arrayNames=arrayNames,
				   arrayInfoColNames=arrayInfo,
				   cdfName="human370v1c",
				   batch=rep("1", nrow(samplesheet)))

}
\dontrun{
	# example for 'krlmm' option
	library(crlmm)
	library(ff)
	# line below is an optional step for krlmm to initialise 16 workers 
	# options("krlmm.cores" = 16)
	# read in raw X and Y intensities output by GenomeStudio's GenCall genotyping module
	XY = readGenCallOutput(c("HumanOmni2-5_4v1_FinalReport_83TUSCAN.csv","HumanOmni2-5_4v1_FinalReport_88CHB-JPT.csv"),
				cdfName="humanomni25quadv1b",
				verbose=TRUE)
	krlmmResult = genotype.Illumina(XY=XY, 
		      			cdfName=ThiscdfName, 
					call.method="krlmm", 
					verbose=TRUE)

	# example for 'krlmm' option with known genotype call for some SNPs and samples
	library(VGAM)
	hapmapCalls = load("hapmapCalls.rda")
	# hapmapCalls should have rownames and colnames corresponding to XY featureNames and sampleNames
	krlmmResult = genotype.Illumina(XY=XY,
					cdfName=ThiscdfName, 
					call.method="krlmm", 
					trueCalls=hapmapCalls, 
					verbose=TRUE)		
}
}
\keyword{classif}

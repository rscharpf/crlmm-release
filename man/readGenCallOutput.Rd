\name{readGenCallOutput}
\alias{readGenCallOutput}

\title{Read X and Y intensities from GenCall output}

\description{
  This function reads the raw X and Y intensities output by GenomeStudio's GenCall genotyping module in preparation for genotyping with \code{crlmm}.
}

\usage{
readGenCallOutput(filenames, path=".", cdfName, colnames=list("SampleID"="Sample ID", "SNPID"="SNP Name", "XRaw"="X Raw", "YRaw"="Y Raw"),
                  type=list("SampleID"="character", "SNPID"="character", "XRaw"="integer", "YRaw"="integer"), verbose=FALSE)
}

\arguments{
  \item{filenames}{'character' string, or a vector of character string specifying the name of the file(s) to read in}
  \item{path}{'character' string specifying the location of file to be read by the function}
  \item{cdfName}{'character' defining the chip annotation (manifest) to use
    ('human370v1c', human550v3b', 'human650v3a', 'human1mv1c',
    'human370quadv3c', 'human610quadv1b', 'human660quadv1a',
    'human1mduov3b', 'humanomni1quadv1b', 'humanomniexpress12v1b', 'humancytosnp12v2p1h')}
  \item{colnames}{list containing elements 'SampleID', 'SNPID', 'XRaw' and 'YRaw', which specify the column names from
     in 'file' that pertain to these variables.  The default should suffice in most situations.}
  \item{type}{list containing data types for the columns to be read in.  The default should be fine in most situations.}
  \item{verbose}{'logical'.  Should processing information be displayed as data is read in?}
}

\details{
This function returns an \code{NChannelSet} containing raw intensity data (X and Y) from GenCall final report file. 
It assumes the GenCall output is formatted to have samples listed one below the other, and that the columns 'X Raw' and 'Y Raw' are available in the file.
The function crlmmillumina() can be run on the output of the \code{readGenCallOutput} function. 
}

\value{
\code{NChannelSet} containing X and Y bead intensities.
}

\references{
  Ritchie ME, Carvalho BS, Hetrick KN, Tavar\'{e} S, Irizarry RA.
  R/Bioconductor software for Illumina's Infinium whole-genome 
  genotyping BeadChips. Bioinformatics. 2009 Oct 1;25(19):2621-3.
}

\author{Cynthia Liu, Matt Ritchie, Zhiyin Dai}

\examples{
#XY = readGenCallOutput(file="Hap650Yv3_Final_Report.txt", cdfName="human650v3a")
#crlmmOut = crlmmIllumina(XY=XY)
}

\keyword{IO}

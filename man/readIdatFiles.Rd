\name{readIdatFiles}
\alias{readIdatFiles}
\alias{readIdatFiles2}

\title{Reads Idat Files from Infinium II Illumina BeadChips}

\description{
  Reads intensity information for each bead type from
  .idat files of Infinium II genotyping BeadChips}

\usage{
readIdatFiles(sampleSheet=NULL, arrayNames=NULL, ids=NULL, path="",
              arrayInfoColNames=list(barcode="SentrixBarcode_A",
                                     position="SentrixPosition_A"),
              highDensity=FALSE, sep="_",
              fileExt=list(green="Grn.idat", red="Red.idat"),
              saveDate=FALSE, verbose=FALSE)
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
  \item{saveDate}{logical.  Should the dates from each .idat be saved
    with sample information?}
  \item{verbose}{logical.  Should processing information be displayed as data is read in?}
}

\details{
The summarised Cy3 (G) and Cy5 (R) intensities (on the orginal scale)
are read in from the .idat files.

Where available, a \code{sampleSheet} data.frame, in the same format
as used by BeadStudio (columns 'Sample\_ID', 'SentrixBarcode\_A' and
'SentrixPosition\_A' are required) which keeps track of sample
information can be specified.

Thanks to Keith Baggerly who provided the code to read in the binary .idat files.
}

\value{
  NChannelSet with intensity data (\code{R}, \code{G}), and indicator 
  for SNPs with 0 beads (\code{zero}) for each bead type.
}

\references{
  Ritchie ME, Carvalho BS, Hetrick KN, Tavar\'{e} S, Irizarry RA.
  R/Bioconductor software for Illumina's Infinium whole-genome 
  genotyping BeadChips. Bioinformatics. 2009 Oct 1;25(19):2621-3.
}

\author{Matt Ritchie}

\examples{

#RG = readIdatFiles()

}
\keyword{IO}

\name{constructAffyCNSet}
\alias{constructAffyCNSet}

\title{ Construct an object of class CNSet from Affymetrix cel files
}
\description{
Construct a container for normalized intensities for Affymetrix cel
files, referred to as a \code{CNSet}
}
\usage{
constructAffyCNSet(filenames, sns, cdfName, batch, verbose = TRUE, genome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filenames}{
    Vector of cel file names.
}
  \item{sns}{
    Sample identifiers. Defaults to \code{basename(filenames)}.
}
  \item{cdfName}{
    Character string indicating annotation package (e.g., "genomewidesnp6Crlmm")
}
  \item{batch}{
    Vector of same length as filenames indicating batch.
}
  \item{verbose}{
    Logical.
}
  \item{genome}{
    Character string indicating UCSC genome build (hg18 or hg19 supported)
}
}

\value{
  An object of class \code{CNSet}
}

\author{
R. Scharpf
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
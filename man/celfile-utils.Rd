\name{validCEL}
\alias{celDates}
\alias{validCEL}
\title{
  Reads cel files and return an error if a file is not read
}
\description{
  Reads cel files and return an error if a file is not read
}
\usage{
validCEL(celfiles)
celDates(celfiles)
}
\arguments{
  \item{celfiles}{
    vector of cel file names to read
}
}
\value{
  Returns a message that cel files were successfully read, or an error
  if there were problems reading the cel files.
}

\author{
R. Scharpf
}

\seealso{ \cite{\link[affyio]{read.celfile.header}},
  \cite{\link{POSIXt}},
  \code{\link{read.celfile}}
}

\examples{
library(oligoClasses)
if(require(hapmapsnp6)){
  path <- system.file("celFiles", package="hapmapsnp6")
  cels <- list.celfiles(path, full.names=TRUE)
  validCEL(cels)
  celDates(cels)
}
}
\keyword{IO}
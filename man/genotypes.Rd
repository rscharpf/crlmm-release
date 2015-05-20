\name{genotypes}
\alias{genotypes}

\title{
  The possible genotypes for an integer copy number.
}
\description{
  The possible genotypes for an integer copy number (0-4).
}
\usage{
genotypes(copyNumber, is.snp=TRUE)
}

\arguments{
  \item{copyNumber}{
    Integer (0-4 allowed).}

  \item{is.snp}{Logical. If TRUE, possible genotypes for a polymorphic
    SNP is returned. If FALSE, only monomorphic genotypes returned.}
}

\value{
  Character vector.
}

\author{
R. Scharpf
}

\examples{

for(i in 0:4) print(genotypes(i))
for(i in 0:4) print(genotypes(i, FALSE))

}

\keyword{manip}
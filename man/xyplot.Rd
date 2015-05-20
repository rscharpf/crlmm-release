\name{xyplot}
\alias{xyplot}
\alias{xyplot,formula,CNSet-method}
\title{Plot prediction regions and normalized intensities.}
\description{
  Plot prediction regions for integer copy number and normalized intensities.
}
\usage{
xyplot(x, data, ...)
}

\arguments{
  \item{x}{
    A \code{formula}.
  }
  \item{data}{
    A \code{CNSet} object.
  }
  \item{\dots}{
    Additional arguments passed to \code{xyplot} function in lattice.
  }
}

\value{
  A \code{trellis} object.
}

\author{
R. Scharpf
}

\seealso{
  \code{\link[lattice]{xyplot}}, \code{\link{ABpanel}}
}

\examples{
library(oligoClasses)
data(cnSetExample2)
table(batch(cnSetExample2))
sample.index <- which(batch(cnSetExample2) == "CUPID")
## A single SNP
pr <- predictionRegion(cnSetExample2[1:4, sample.index], copyNumber=0:4)
gt <- calls(cnSetExample2[1:4, sample.index])
lim <- c(6,13)
xyplot(B~A|snpid, data=cnSetExample2[1:4, sample.index],
       predictRegion=pr,
       panel=ABpanel,
       pch=21,
       fill=c("red", "blue", "green3")[gt],
       xlim=lim, ylim=lim)

## multiple SNPs, prediction regions for 3 batches
\dontrun{
	tab <- table(batch(cnSetExample2))
	bns <- names(tab)[tab > 50]
	sample.index <- which(batch(cnSetExample2) %in% bns[1:3])
	pr <- predictionRegion(cnSetExample2[1:10, sample.index], copyNumber=0:4)
	gt <- as.integer(calls(cnSetExample2[1:10, sample.index]))
	xyplot(B~A|snpid, data=cnSetExample2[1:10, sample.index],
	       predictRegion=pr,
	       panel=ABpanel,
	       pch=21,
	       fill=c("red", "blue", "green3")[gt],
	       xlim=c(6,12), ylim=c(6,12))

	## nonpolymorphic markers
	data(cnSetExample2)
	tab <- table(batch(cnSetExample2))
	bns <- names(tab)[tab > 50]
	sample.index <- which(batch(cnSetExample2)%in%bns[1:3])
	np.index <- which(!isSnp(cnSetExample2))[1:10]
	taus <- tau2(cnSetExample)[np.index, , , ]
	pr <- predictionRegion(cnSetExample2[np.index, sample.index],
			       copyNumber=0:4)
	pp <- posteriorProbability(cnSetExample2[np.index, sample.index],
				   predictRegion=pr,
				   copyNumber=0:4)
}
}
\keyword{dplot}
\keyword{hplot}

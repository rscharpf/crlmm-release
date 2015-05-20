\name{ABpanel}
\alias{ABpanel}
\title{A panel function for plotting prediction regions and
  log-normalized intensities
}
\description{
  A panel function for plotting prediction regions and
  log-normalized intensities
}
\usage{
ABpanel(x, y, predictRegion, copyNumber = 0:4, fill, ..., subscripts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    log-normalized intensities for the A or B allele
}
\item{y}{
      log-normalized intensities for the A or B allele
}
  \item{predictRegion}{
    A \code{list}. See \code{predictionRegion}.
}
  \item{copyNumber}{
    Integer vector.  Indicates which prediction regions are drawn.
}
  \item{fill}{
    Character or integer vector for coloring the points.  Only valid for
    certain point symbols.  See \code{points}.
}
  \item{\dots}{
    Additional arguments to \code{panel.xyplot} and \code{lpolygon}.
}
  \item{subscripts}{
    See \code{xyplot} in the \pkg{lattice} package.
}
}

\value{
  Not applicable
}

\note{

  \code{ABpanel} can be passed as the argument to panel in the
  \code{xyplot} method for \code{CNSet} objects.  See the examples in
  \code{xyplot}.

}

\author{
R. Scharpf
}

\seealso{
  \code{\link[lattice]{xyplot}},   \code{\link[lattice]{panel.xyplot}}
  \code{\link[lattice]{lpolygon}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{aplot}
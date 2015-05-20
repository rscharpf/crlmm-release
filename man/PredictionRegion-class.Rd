\name{PredictionRegion-class}
\Rdversion{1.1}
\docType{class}
\alias{PredictionRegion-class}
\alias{[,PredictionRegion,ANY,ANY,ANY-method}

\title{Class \code{"PredictionRegion"}}
\description{
  A container for bivariate normal prediction regions for SNP data and
  univarite prediction regions for nonpolymorphic markers.
}
\section{Objects from the Class}{
  Objects from the class are created from the \code{predictionRegion} function.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
Class \code{"\linkS4class{list_or_ffdf}"}, by class "list", distance 2.
Class \code{vectorORfactor}, by class "list", distance 3.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "PredictionRegion")}: ... } Prediction
    regions can be subset by markers.
	 }
}

\author{
R. Scharpf
}


\seealso{
\code{\link{predictionRegion}}
}
\examples{
showClass("PredictionRegion")
}
\keyword{classes}

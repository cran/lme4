\name{summary.lmer-class}
\docType{class}
\alias{summary.lmer-class}
\alias{show,summary.lmer-method}

\title{Class "summary.lmer"}
\description{A summary of an \code{\link{lmer}} fitted model object.}
\section{Objects from the Class}{
  Objects are created by applying the \code{summary} generic to an
  \code{\link{lmer}} object.
}
\section{Slots}{
  \describe{
    \item{\code{useScale}:}{logical indicator of the use of the scale
      factor when calculating variance-covariance objects.}
    \item{\code{showCorrelation}:}{Logical indicator of whether to
      display the correlation of the fixed-effects estimates when
      showing the object.}
    All the other slots are inherited from the
    \code{\link{lmer-class}{lmer}} class which is contained by
    this class.
  }
}
\section{Extends}{
Class \code{"lmer"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "summary.lmer")}:}
  }
}
\seealso{
  \code{\link[Matrix]{lmer-class}}
}
%\examples{}
\keyword{classes}

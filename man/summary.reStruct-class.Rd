\name{summary.reStruct-class}
\docType{class}
\alias{summary.reStruct-class}
\title{Class "summary.reStruct"}
\description{A summary of an \code{\link[lme4:reStruct-class]{reStruct}}
  object.  The summary contains information derived from the original
  object.} 
\section{Objects from the Class}{
  Objects from the class are usually constructed by the
  \code{\link[base]{summary}} generic function applied to an
  \code{\link[lme4:reStruct-class]{reStruct}} object.
}
\section{Slots}{
  \describe{
    \item{\code{coefficients}:}{A matrix containing the fixed-effects
      coefficients, their standard errors and the t statistics.}
    \item{\code{sigma}:}{The estimated standard deviation of the lowest
      level noise term in the model.}
    \item{\code{logLik}:}{The log-likelihood at the estimated parameter
      values.}
  }
}

\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "summary.reStruct")}: print the object}
  }
}
\seealso{\code{\link{reStruct-class}}, \code{\link[base]{summary}}}
%\examples{}
\keyword{classes}

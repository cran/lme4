\name{summary.lme-class}
\docType{class}
\alias{summary.lme-class}
\alias{summary.ssclme-class}
\alias{coef,summary.lme-method}
\alias{coef,summary.ssclme-method}
\alias{show,summary.lme-method}
\alias{show,summary.ssclme-method}
\title{Class "summary.lme"}
\description{A summary of an \code{"lme"} fitted model object.}
\section{Objects from the Class}{
  Objects are created by applying the \code{summary} generic to an
  \code{"lme"} object.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{A copy of the function call that created the
      \code{"lme"} object.}
    \item{\code{logLik}:}{A \code{"logLik"} object for the log-likelihood.}
    \item{\code{re}:}{A \code{"summary.ssclme"} object summarizing the
      model representation.} 
    \item{\code{residuals}:}{The residuals as a \code{"numeric"} vector.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "summary.lme")}: Display the object}
  }
}
\seealso{
  \code{\link{lme-class}}, \code{\link[Matrix]{ssclme-class}}
}
\keyword{classes}

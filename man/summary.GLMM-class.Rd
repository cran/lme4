\name{summary.GLMM-class}
\docType{class}
\alias{summary.GLMM-class}
\alias{show,summary.GLMM-method}
\title{Class "summary.GLMM"}
\description{A summary of an \code{"GLMM"} fitted model object.}
\section{Objects from the Class}{
  Objects are created by applying the \code{summary} generic to an
  \code{"GLMM"} object.
}
\section{Slots}{
  \describe{
    \item{\code{family}:}{The family of functions
      defining the generalized linear model.  See
      \code{\link[stats]{family}} for details.}
    \item{\code{call}:}{A copy of the function call that created the
      \code{"GLMM"} object.}
    \item{\code{logLik}:}{A \code{"logLik"} object for the log-likelihood.}
    \item{\code{re}:}{A \code{"summary.ssclme"} object summarizing the
      model representation.} 
    \item{\code{residuals}:}{The residuals as a \code{"numeric"} vector.}
  }
}
\section{Extends}{
Class \code{"summary.lme"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "summary.GLMM")}: Display the object.}
  }
}
\seealso{
  \code{\link{GLMM}}, \code{\link{GLMM-class}},
  \code{\link[Matrix]{ssclme-class}}
}
\keyword{classes}

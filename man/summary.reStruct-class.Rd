\name{summary.reStruct-class}
\docType{class}
\alias{summary.reStruct-class}
\alias{show,summary.reStruct-method}
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
    \item{\code{fixed}:}{\code{"formula"}: the formula defining the
      fixed-effects.}
    \item{\code{coefficients}:}{A matrix containing the fixed-effects
      coefficients, their standard errors and the t statistics.}
    \item{\code{scale}:}{\code{"numeric"}.  The estimate of the scale
      parameter in the distribution of the per-observation random
      variation, if it is estimated.  In the case of a \code{"lme"}
      model this is the estimated standard deviation of the lowest
      level noise term in the model.}
    \item{\code{denomDF}:}{\code{"integer"}: vector of denominator
      degrees of freedom.}
    \item{\code{REML}:}{\code{"logical"}: was the REML criterion used to
      estimate the parameters?}
    \item{\code{random}:}{\code{"lmeLevelList"}: as in the
      \code{\link[lme4:reStruct-class]{reStruct}} object from which this
      was derived.}
    \item{\code{ngrps}:}{\code{"integer"}: number of groups at each
      level of random effects.}
    \item{\code{nobs}:}{\code{"integer"}: total number of observations.}
    \item{\code{corFixed}:}{\code{"corrmatrix"}: Correlation of the
      estimates of the fixed-effects parameters.}
    \item{\code{reSumry}:}{\code{"list"}: a summary of the random
      effects parameter estimates}
    \item{\code{useScale}:}{\code{"logical"}: should the \code{"scale"}
      slot be used when calculating the variances and covariances of the
      random effects?}
    \item{\code{showCorrelation}:}{\code{"logical"}: should the
      \code{"corFixed"} slot be shown in the summary.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(x = "summary.reStruct")}: show the object}
  }
}
\seealso{\code{\link{reStruct-class}}, \code{\link[base]{summary}}}
%\examples{}
\keyword{classes}

\name{GLMM-class}
\docType{class}
\alias{GLMM-class}
\alias{coef,GLMM-method}
\alias{deviance,GLMM-method}
\alias{fixef,GLMM-method}
\alias{logLik,GLMM-method}
\alias{ranef,GLMM-method}
\alias{show,GLMM-method}
\alias{summary,GLMM-method}

\title{Class "GLMM"}
\description{A fitted generalized linear mixed model.}
\section{Objects from the Class}{
    Objects are usually created by calls to the constructor function
  \code{GLMM}.  They also can be created by calls of the form
  \code{new("GLMM", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{family}:}{The family of functions
      defining the generalized linear model.  See
      \code{\link[stats]{family}} for details.}
    \item{\code{logLik}:}{the \code{"numeric"} value of the log-likelihood.}
    \item{\code{fixef}:}{A \code{"numeric"} vector of fixed-effects
      parameters.}
    \item{\code{call}:}{A copy of the function call that created the object.}
    \item{\code{facs}:}{A \code{list} of (possibly reordered) grouping
      factors associated with the random effects.}
    \item{\code{x}:}{If the optional argument \code{x} to
      \code{lme} is \code{TRUE}, a \code{list} of model matrices
      associated with the random effects, and the fixed effects with the
      response appended. Otherwise, an empty list.}
    \item{\code{model}:}{The model frame (of class \code{"data.frame"})
      for the model or, if the optional argument
      \code{model} to \code{lme} is \code{FALSE}, an empty frame.}
    \item{\code{REML}:}{A \code{"logical"} indicator of the model having
      been fit according to the REML criterion.}
    \item{\code{rep}:}{An \code{"ssclme"} object representing the fitted
      model.}
    \item{\code{fitted}:}{A \code{"numeric"} vector of fitted values.}
    \item{\code{residuals}:}{A \code{"numeric"} vector of raw residuals.}
  }
}
\section{Extends}{
Class \code{"lme"}, directly.
}
\section{Methods}{
  \describe{
    \item{VarCorr}{\code{signature(x = "GLMM")}: Extract the variances,
      standard deviations, and correlations of the random effects.}
    \item{anova}{\code{signature(object = "GLMM")}: Perform an analysis
      of variance.}
    \item{coef}{\code{signature(object = "GLMM")}: Extract the parameters
      that determine the relative precision matrices. The optional
      argument \code{unconst} determines if the constrained or
      unconstrained parameterization is used.}
    \item{deviance}{\code{signature(object = "GLMM")}: Extract the
      deviance as a numeric scalar.  The optional argument \code{REML}
      determines if the REML or ML criterion is used.}
    \item{fitted}{\code{signature(object = "GLMM")}: Extract the fitted
      values as a numeric vector.}
    \item{fixef}{\code{signature(object = "GLMM")}: Extract the fixed
      effects coefficients as a named numeric vector }
    \item{formula}{\code{signature(x = "GLMM")}: Extract the formula of
      the response and the fixed effects.}
    \item{logLik}{\code{signature(object = "GLMM")}: Extract the
      log-likelihood.}
    \item{plot}{\code{signature(x = "GLMM")}:}
    \item{ranef}{\code{signature(object = "GLMM")}: Extact the random
      effects as a named list of numeric matrices.}
    \item{residuals}{\code{signature(object = "GLMM")}: Extract the
      residuals as a numeric vector.}
    \item{show}{\code{signature(object = "GLMM")}: Print a concise
      description of the object.}
    \item{summary}{\code{signature(object = "GLMM")}: Create a summary
      object of class \code{"summary.lme"}.}
    \item{update}{\code{signature(object = "GLMM")}: Create an updated
      fitted model.}
    \item{vcov}{\code{signature(object = "GLMM")}: Extract the variances
      and covariances of the fixed-effects parameter estimates.}
  }
}
\details{Many methods for the \code{"GLMM"} class simply recall the
  generic on \code{as(object, "lme")} which, in turn, often calls the
    generic on the \code{"rep"} component of the first argument.  Thus the
    methods actually are applied to an object of the \code{"lme"} or
    \code{"ssclme"} class.}
\seealso{
  \code{\link{GLMM}}, \code{\link{lme-class}},
  \code{\link[Matrix]{ssclme-class}}
}
\keyword{classes}

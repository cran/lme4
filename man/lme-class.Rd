\name{lme-class}
\docType{class}
\alias{lme-class}
\alias{VarCorr,lme-method}
\alias{VarCorr,ssclme-method}
\alias{anova,lme-method}
\alias{anova,ssclme-method}
\alias{coef,lme-method}
\alias{deviance,lme-method}
\alias{deviance,ssclme-method}
\alias{fitted,lme-method}
\alias{formula,lme-method}
\alias{logLik,lme-method}
\alias{plot,lme-method}
\alias{ranef,lme-method}
\alias{ranef,ssclme-method}
\alias{residuals,lme-method}
\alias{show,lme-method}
\alias{show,ssclme-method}
\alias{summary,lme-method}
\alias{summary,ssclme-method}
\alias{update,lme-method}
\alias{vcov,lme-method}
\alias{vcov,ssclme-method}
\title{Class "lme"}
\description{A fitted linear mixed-effects model.}
\section{Objects from the Class}{
  Objects are usually created by calls to the constructor function
  \code{lme}.  They also can be created by calls of the form
  \code{new("lme", ...)}.
}
\section{Slots}{
  \describe{
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
\section{Methods}{
  \describe{
    \item{VarCorr}{\code{signature(x = "lme")}: Extract the variances,
      standard deviations, and correlations of the random effects.}
    \item{anova}{\code{signature(object = "lme")}: Perform an analysis
      of variance.}
    \item{coef}{\code{signature(object = "lme")}: Extract the parameters
      that determine the relative precision matrices. The optional
      argument \code{unconst} determines if the constrained or
      unconstrained parameterization is used.}
    \item{deviance}{\code{signature(object = "lme")}: Extract the
      deviance as a numeric scalar.  The optional argument \code{REML}
      determines if the REML or ML criterion is used.}
    \item{fitted}{\code{signature(object = "lme")}: Extract the fitted
      values as a numeric vector.}
    \item{fixef}{\code{signature(object = "lme")}: Extract the fixed
      effects coefficients as a named numeric vector }
    \item{formula}{\code{signature(x = "lme")}: Extract the formula of
      the response and the fixed effects.}
    \item{logLik}{\code{signature(object = "lme")}: Extract the
      log-likelihood corresponding to the ML or REML criterion.}
    \item{plot}{\code{signature(x = "lme")}:}
    \item{ranef}{\code{signature(object = "lme")}: Extact the random
      effects as a named list of numeric matrices.}
    \item{residuals}{\code{signature(object = "lme")}: Extract the
      residuals as a numeric vector.}
    \item{show}{\code{signature(object = "lme")}: Print a concise
      description of the object.}
    \item{summary}{\code{signature(object = "lme")}: Create a summary
      object of class \code{"summary.lme"}.}
    \item{update}{\code{signature(object = "lme")}: Create an updated
      fitted model.}
    \item{vcov}{\code{signature(object = "lme")}: Extract the variances
      and covariances of the fixed-effects parameter estimates.}
  }
}
\details{Many methods for the \code{"lme"} class simply recall the
  generic on the \code{"rep"} component of the first argument.  Thus the
  methods actually are applied to an object of the \code{"ssclme"} class.}
\seealso{
  \code{\link{lme}}, \code{\link[Matrix]{ssclme-class}}
}
\keyword{classes}

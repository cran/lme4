\name{lme-class}
\docType{class}
\alias{lme-class}
\alias{VarCorr,lme-method}
\alias{anova,lme-method}
\alias{coef,lme-method}
\alias{fitted,lme-method}
\alias{fixef,lme-method}
\alias{formula,lme-method}
\alias{getResponse,lme-method}
\alias{intervals,lme-method}
\alias{logLik,lme-method}
\alias{plot,lme-method}
\alias{ranef,lme-method}
\alias{residuals,lme-method}
\alias{show,lme-method}
\alias{summary,lme-method}
\alias{update,lme-method}
\title{Class "lme"}
\description{A linear mixed-effects model.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("lme", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{reStruct}:}{Object of class \code{"reStruct"},
      the random-effects structure for the model.}
    \item{\code{frame}:}{Object of class \code{"data.frame"},
      the model.frame used to fit the model.}
    \item{\code{na.action}:}{Object of class \code{"ANY"},
      the \code{na.action} argument used when creating the
      \code{model.frame}.}
    \item{\code{fitted}:}{Object of class \code{"numeric"}, the current
      fitted values.}
    \item{\code{call}:}{Object of class \code{"call"},
      the function call used to create the object.}
  }
}
\section{Extends}{
Class \code{"lme"}, directly.
}
\section{Methods}{
  \describe{
    \item{VarCorr}{\code{signature(x = "lme")}: extract variance and correlation components.}
    \item{anova}{\code{signature(object = "lme")}: calculate analysis of
      variance tables for linear mixed effects models.}
    \item{coef}{\code{signature(object = "lme")}: extract the
      random effects parameters for the model.}
    \item{fitted}{\code{signature(object = "lme")}: extract the fitted values.}
    \item{fixef}{\code{signature(object = "lme")}: extract the fixed
      effects parameters.}
    \item{formula}{\code{signature(x = "lme")}: extract the fixed
      effects formula.}
    \item{getResponse}{\code{signature(object = "lme")}: extract the response.}
    \item{intervals}{\code{signature(object = "lme")}: calculate
      confidence intervals.}
    \item{logLik}{\code{signature(object = "lme")}: calculate the
      log-likelihood or log-restricted-likelihood for the model}
    \item{plot}{\code{signature(x = "lme")}: plot the fitted model.}
    \item{ranef}{\code{signature(object = "lme")}: calculate the
      random-effects coefficients}
    \item{residuals}{\code{signature(object = "lme")}: calculate the
      residuals from the model fit.}
    \item{show}{\code{signature(object = "lme")}: show the object.}
    \item{summary}{\code{signature(object = "lme")}: summarize the object.}
    \item{update}{\code{signature(object = "lme")}: update the model fit.}
  }
}
\author{Saikat DebRoy \email{saikat@stat.wisc.edu}
  and Douglas Bates \email{bates@stat.wisc.edu}}
\examples{
library(lme4)
data(bdf)
# This returns an object of class lme
fm <- lme(langPOST ~ IQ.ver.cen + avg.IQ.ver.cen, data = bdf,
          random = ~ IQ.ver.cen | schoolNR)
fm
}
\keyword{classes}

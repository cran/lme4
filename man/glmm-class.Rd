\name{glmm-class}
\docType{class}
\alias{glmm-class}
\alias{VarCorr,glmm-method}
\alias{fixef<-,glmm,numeric-method}
\alias{getResponse,glmm-method}
\alias{logLik,glmm-method}
\alias{show,glmm-method}
\alias{summary,glmm-method}
\title{Class "glmm"}
\description{A generalized linear mixed-effects model.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("glmm", ...)}
  but more commonly they are created by calls to \code{\link{GLMM}}.
}
\section{Slots}{
  \describe{
    \item{\code{family}:}{Object of class \code{"family"}, specifying
      the error distribution and the link function used in the model.}
    \item{\code{origy}:}{Object of class \code{"numeric"}, the original
      response.}
    \item{\code{n}:}{Object of class \code{"numeric"}, if the family is
      \code{binomial}, the number of observations for each data point.}
    \item{\code{prior.weights}:}{Object of class \code{"numeric"},
      weights used when creating the model.frame.}
    \item{\code{init.weights}:}{Object of class \code{"numeric"},
      weights from a fixed effects generalized linear model for the data.}
    \item{\code{init.y}:}{Object of class \code{"numeric"}, the linear
      predictor from a fixed effects generalized linear model for the data.}
    \item{\code{method}:}{Object of class \code{"character"}, the method
      used to fit the generalized linear mixed model.}
    \item{\code{reStruct}:}{Object of class \code{"reStruct"}, from class
	\code{"lme"}, the random-effects structure for the model.}
    \item{\code{frame}:}{Object of class \code{"data.frame"}, from class
	\code{"lme"}, the model.frame used to fit the model.}
    \item{\code{na.action}:}{Object of class \code{"ANY"}, from class
	\code{"lme"}, the \code{na.action} argument used when creating the
      \code{model.frame}.}
    \item{\code{fitted}:}{Object of class \code{"numeric"}, from class
      \code{"lme"}, the fitted values in the linear predictor scale.}
    \item{\code{call}:}{Object of class \code{"call"}, from class
      \code{"lme"}, the function call used to create the object.}
  }
}
\section{Extends}{
Class \code{"lme"}, directly.
}
\section{Methods}{
  \describe{
    \item{LMEoptimize<-}{\code{signature(x = "glmm", value = "list")}:
      optimize the PQL approximation to the log-likelihood.}
    \item{VarCorr}{\code{signature(x = "glmm")}: Extract variance and correlation components.}
    \item{fixef<-}{\code{signature(object = "glmm", value = "numeric")}:
      assign the fixed effects (used for \code{method = "Laplace"}).}
    \item{getResponse}{\code{signature(object = "glmm")}: extract the response.}
    \item{logLik}{\code{signature(object = "glmm")}: extract the
      (approximate) log-likelihood.}
    \item{show}{\code{signature(object = "glmm")}: show the object.}
    \item{summary}{\code{signature(object = "glmm")}: summarize the object.}
  }
}
\author{Saikat DebRoy \email{saikat@stat.wisc.edu}
  and Douglas Bates \email{bates@stat.wisc.edu}}
\examples{
library(lme4)
data(guImmun)
# This returns an object of class glmm
fm = GLMM(immun ~ kid2p + mom25p + ord + ethn +
                  momEd + husEd + momWork + rural + pcInd81,
          data = guImmun, family = binomial,
          random = ~1|comm/mom)
fm
}
\keyword{classes}

\name{glmm-class}
\docType{class}
\alias{glmm-class}
\title{Class "glmm"}
\description{A generalized linear mixed-effects model.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("glmm", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{family}:}{Object of class \code{"family"}, specifying
      the error distribution and the link function used in the model.}
    \item{\code{origy}:}{Object of class \code{"numeric"}, the original
      response.}
    \item{\code{n}:}{Object of class \code{"numeric"}, if the family is
      \code{binomial}, the number of observations for each data point.}
    \item{\code{offset}:}{Object of class \code{"numeric"}, offset when
      specifying the model.}
    \item{\code{prior.weights}:}{Object of class \code{"numeric"},
      weights used when creating the model.frame.}
    \item{\code{init.weights}:}{Object of class \code{"numeric"},
      weights from a fixed effects generalized linear model for the data.}
    \item{\code{init.y}:}{Object of class \code{"numeric"}, the linear
      predictor from a fixed effects generalized linear model for the data.}
    \item{\code{method}:}{Object of class \code{"character"}, the method
      used to fit the generalized linear mixed model.}
    \item{\code{reStruct}:}{Object of class \code{"reStruct", from class
	"lme"}, the random-effects structure for the model.}
    \item{\code{frame}:}{Object of class \code{"data.frame", from class
	"lme"}, the model.frame used to fit the model.}
    \item{\code{na.action}:}{Object of class \code{"ANY", from class
	"lme"}, the \code{na.action} argument used when creating the
      \code{model.frame}.}
    \item{\code{call}:}{Object of class \code{"call", from class "lme"},
      the function call used to create the object.}
  }
}
\section{Extends}{
Class \code{"lme"}, directly.
}
\section{Methods}{
  \describe{
    \item{VarCorr}{\code{signature(x = "glmm")}: Extract variance and correlation components.}
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

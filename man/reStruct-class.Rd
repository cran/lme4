\name{reStruct-class}
\docType{class}
\alias{reStruct-class}
\alias{LMEgradient,reStruct,missing,missing-method}
\alias{VarCorr,reStruct-method}
\alias{coef,reStruct-method}
\alias{fitted,reStruct-method}
\alias{fixef<-,reStruct,numeric-method}
\alias{getFixDF,reStruct-method}
\alias{getGroups,reStruct,missing,ANY,missing,missing-method}
\alias{getResponse,reStruct-method}
\alias{logLik,reStruct-method}
\alias{model.matrix,reStruct-method}
\alias{model.matrix<-,reStruct-method}
\alias{response<-,reStruct,numeric-method}
\alias{summary,reStruct-method}
\alias{weighted<-,reStruct,numeric-method}
\title{Class "reStruct"}
\description{The random-effects model structure in a linear
  mixed-effects model or a generalized linear mixed-effects model or a
  nonlinear mixed-effects model.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("reStruct", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{random}:}{A \code{"list"} of \code{lmeLevel} objects
      giving the levels of random effects in the model, the
      fixed-effects, and the response or the working residual.}
    \item{\code{fixed}:}{The \code{"formula"} of the response and
      fixed-effects terms.}
    \item{\code{offset}:}{\code{"numeric"}: the offset in the model, if present.}
    \item{\code{dirtyDecomposed}:}{\code{"logical"}: if \code{TRUE} the
      decomposed matrix must be recalculated.}
    \item{\code{useWeighted}:}{\code{"logical"}: if \code{TRUE}
      calculate and use the weighted model matrix to form the
      decompositions; otherwise use the original model matrix.}
    \item{\code{dirtyStored}:}{\code{"logical"}: if \code{TRUE} the
      matrix of stored decompositions must be recalculated.}
    \item{\code{dirtyBbetas}:}{\code{"logical"}: if \code{TRUE} the
      BLUP's and the conditional estimates of the fixed-effects must be
      recalculated.}
    \item{\code{logLik}:}{\code{"numeric"}: the log-likelihood at the
      current parameter values.}
    \item{\code{REML}:}{\code{"logical"}: if \code{TRUE} the parameters
      will be estimated according to the REML criterion}
    \item{\code{analyticHessian}:}{\code{"logical"}: if \code{TRUE} the
      \code{hessianArray} slot of the lmeLevel objects are calculated
      during optimization.}
    \item{\code{reverseOrder}:}{\code{"integer"}: the permutation of the
      rows that provides the original order.}
    \item{\code{origOrder}:}{\code{"integer"}: the inverse of the
      \code{reverseOrder} permutation.}
    \item{\code{original}:}{\code{"matrix"}: the original model matrix,
      including the original response in the last column.}
    \item{\code{weighted}:}{\code{"matrix"}: the model matrix and
      responses after applying weights.}
    \item{\code{stored}:}{\code{"matrix"}: a model matrix of
      intermediate decomposition results needed for evaluating the BLUPs
      and the EM or ECME iterations.}
    \item{\code{decomposed}:}{\code{"matrix"}: the model matrix after
      predecomposition.  This generally has many few rows than
      \code{original}.}
    \item{\code{bbetas}:}{\code{"numeric"}: the BLUPs and the conditional
       estimates of the fixed-effects parameters at the current values
       of the relative precision matrices.}
    \item{\code{dontCopy}:}{\code{"logical"}: if \code{TRUE} it indicates
      that this object has just been created and is only assigned to one
      name.  In these circumstances changes are made directly on the
      object without copying. This is dangerous.  You probably don't
      want to modify this setting.}
    \item{\code{assign.X}:}{Object of class \code{"ANY"}: the
      \code{assign} attribute from the model matrix for the fixed effects.}
  }
}
\section{Methods}{
  \describe{
    \item{EMsteps<-}{\code{signature(x = "reStruct", value = "list")}:
      perform the EM iterations}
    \item{LMEgradient}{\code{signature(x = "reStruct", A = "missing",
	nlev = "missing")}: evaluate the gradient of the linear
      mixed-effects profiled log-likelihood.}
    \item{LMEoptimize<-}{\code{signature(x = "reStruct", value =
	"list")}: Optimize with \code{optim} or \code{nlm}.}
    \item{VarCorr}{\code{signature(x = "reStruct")}: Extract the
      variances and covariances of the random effects.}
    \item{coef}{\code{signature(object = "reStruct")}: return the
      current parameters in the object.}
    \item{coef<-}{\code{signature(object = "reStruct", value =
	"numeric")}: assign the parameters to the object.}
    \item{fitted}{\code{signature(object = "reStruct")}: extract the
      fitted values.}
    \item{fixef}{\code{signature(object = "reStruct")}: extract the
      fixed effects.}
    \item{fixef<-}{\code{signature(object = "reStruct", value =
	"numeric")}: assign the fixed effects.}
    \item{getFixDF}{\code{signature(object = "reStruct")}: extract the
      denominator degrees of freedom for the fixed-effects terms.}
    \item{getGroups}{\code{signature(object = "reStruct", form =
	"missing", level = "ANY", data = "missing", sep = "missing")}:
      return the grouping factors for the random effects}
    \item{getResponse}{\code{signature(object = "reStruct")}: extract
      the response.}
    \item{logLik}{\code{signature(object = "reStruct")}: return the
      log-likelihood at the current parameter values. }
    \item{model.matrix}{\code{signature(object = "reStruct")}: return the model
      matrix for the object.}
    \item{model.matrix<-}{\code{signature(x = "reStruct", value =
	"matrix")}: assign the model matrix to the object.}
    \item{ranef}{\code{signature(object = "reStruct")}: extract the
      random effects predictors.}
    \item{response<-}{\code{signature(x = "reStruct", value =
	"numeric")}: set the working response (used in GLMM models).}
    \item{summary}{\code{signature(object = "reStruct")}: summarize the object.}
    \item{weighted<-}{\code{signature(x = "reStruct", value =
	"matrix")}: update the weights.}
  }
}
%\references{}
%\author{Saikat DebRoy \email{saikat@stat.wisc.edu} and Douglas Bates \email{bates@stat.wisc.edu}}
%\note{}
%\seealso{}
\examples{}
\keyword{classes}

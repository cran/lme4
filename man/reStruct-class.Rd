\name{reStruct-class}
\docType{class}
\alias{reStruct-class}
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
    \item{\code{dirty}:}{\code{"logical"}: if \code{TRUE} the
      decomposition must be updated. }
    \item{\code{logLik}:}{\code{"numeric"}: the log-likelihood at the
      current parameter values.}
    \item{\code{frame}:}{A model frame in which the \code{fixed} formula
      is evaluated.}
    \item{\code{REML}:}{\code{"logical"}: if \code{TRUE} the parameters
      will be estimated according to the REML criterion}
    \item{\code{EMiter}:}{\code{"integer"}: the number of EM iterations
      to be performed.}
    \item{\code{groups}:}{\code{"data.frame"}: the grouping factors for
      the random-effects terms.}
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
    \item{\code{reverseOrder}:}{\code{"integer"}: the permutation of the
      rows that provides the original order.}
    \item{\code{bbetas}:}{\code{"numeric"}: the BLUPs and the conditional
       estimates of the fixed-effects parameters at the current values
       of the relative precision matrices.}
  }
}
\section{Methods}{
  \describe{
    \item{EMsteps<-}{\code{signature(x = "reStruct", value = "list")}:
      perform the EM iterations}
    \item{LMEgradient}{\code{signature(x = "reStruct", A = "missing",
	nlev = "missing")}: evaluate the gradient of the linear
      mixed-effects log-likelihood.}
    \item{LMEoptimize<-}{\code{signature(x = "reStruct", value =
	"list")}: Optimize with \code{optim} or \code{nlm}.}
    \item{coef}{\code{signature(object = "reStruct")}: return the
      current parameters in the object.}
    \item{coef<-}{\code{signature(object = "reStruct", value =
	"numeric")}: assign the parameters to the object.}
    \item{groups}{\code{signature(x = "reStruct")}: return the grouping
      factors for the random effects} 
    \item{logLik}{\code{signature(object = "reStruct")}: return the
      log-likelihood at the current parameter values. }
    \item{model.matrix}{\code{signature(object = "reStruct")}: return the model
      matrix for the object.}
    \item{model.matrix<-}{\code{signature(x = "reStruct", value =
	"matrix")}: assign the model matrix to the object.}
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

\name{lmeLevel-class}
\docType{class}
\alias{lmeLevel-class}
\alias{coef,lmeLevel-method}
\title{Class "lmeLevel"}
\description{An object representing a single level of random effects in
  a linear mixed-effects model.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("lmeLevel",
    ...)} or by calls to the generic constructor \code{lmeLevel}.
}
\section{Slots}{
  \describe{
    \item{\code{precision}:}{A precision matrix of class \code{"pdMat"}}
    \item{\code{groups}:}{The grouping factor (of class \code{"factor"}) for the level.}
    \item{\code{columns}:}{An \code{"integer"} vector of the columns in
      the model matrix that correspond to this level.}
    \item{\code{parsInd}:}{An \code{"integer"} vector giving the
      subscripts in the overall parameter vector that apply to this level.}
    \item{\code{originalRows}:}{A list of integer vectors of row
      numbers in the original model matrix that correspond to each of the
      levels of the grouping factor.}
    \item{\code{decomposedRows}:}{A list of integer vectors of row
      numbers in the decomposed model matrix that correspond to each of the
      levels of the grouping factor.}
    \item{\code{storedRows}:}{A list of integer vectors of row
      numbers in the stored data matrix that correspond to each of the
      levels of the grouping factor.}
    \item{\code{nrow}:}{\code{integer}: the number of rows in the
      precision matrix.}
    \item{\code{updateFactor}:}{Object of class \code{"matrix"} - the
      matrix \code{A} used to evaluate the gradient of the
      log-likelihood and to formulate the ECME step.}
    \item{\code{hessianArray}:}{Object of class \code{"array"} - the
      array \code{H} used to evaluate the Hessian of the
      log-likelihood.}
    \item{\code{nlev}:}{\code{integer}: the number of levels in the
      grouping factor}
  }
}
\section{Methods}{
  \describe{
    \item{EMupdate<-}{\code{signature(x = "lmeLevel", nlev = "missing",
	value = "matrix")}: update the level within the EM or ECME algorithm}
    \item{LMEgradient}{\code{signature(x = "lmeLevel", A = "matrix",
	nlev = "missing")}: calculate the gradient of the log-likelihood
      with respect to the parameters from the decomposition of the scaled
      BLUPs and their variances.}
    \item{coef}{\code{signature(object = "lmeLevel")}: extract the
      parameters (from the precision matrix)}
    \item{coef<-}{\code{signature(object = "lmeLevel", value =
	"numeric")}: set the parameters.}
  }
}
\note{In \code{reStruct} objects the random-effects terms and the
  fixed-effects terms and the response are all represented by
  \code{lmeLevel} objects.  The \code{precision} slot in the
  fixed-effects and response objects is an uninitialized \code{pdMat}
  objects and the \code{nrow} slot is zero.}
%\references{}
\author{Saikat DebRoy \email{saikat@stat.wisc.edu} and Douglas Bates
  \email{bates@stat.wisc.edu}}
\seealso{\code{\link{lmeLevel}}}
%\examples{}
\keyword{classes}

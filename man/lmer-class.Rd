\name{lmer-class}
\docType{class}
\alias{lmer-class}
\alias{VarCorr,lmer-method}
\alias{anova,lmer-method}
\alias{coerce,lmer,dsTMatrix-method}
\alias{coerce,lmer,dtTMatrix-method}
\alias{coef,lmer-method}
\alias{confint,lmer-method}
\alias{deviance,lmer-method}
\alias{formula,lmer-method}
\alias{fitted,lmer-method}
\alias{logLik,lmer-method}
\alias{residuals,lmer-method}
\alias{resid,lmer-method}
\alias{show,lmer-method}
\alias{summary,lmer-method}
\alias{update,lmer-method}
\alias{vcov,lmer-method}
\alias{plot,lmer.ranef-method}
\alias{plot,lmer.coef-method}
\alias{terms,lmer-method}
\alias{with,lmer-method}
\title{lmer model representation}
\description{A representation of a linear mixed-effects model using a
  blocked, sparse, symmetric, column-oriented matrix and dense matrices.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("lmer", ...)}
or, more commonly, via the \code{lmer} function.}
}
\section{Slots}{
  \describe{
     \item{\code{flist}:}{The list of (possibly permuted) grouping
       factors for the random effects.}
    \item{\code{perm}:}{A list of permutation vectors applied to the
      grouping factors.}
    \item{\code{Parent}:}{The extended parent array for the diagonal
      elements of \code{L}.}
    \item{\code{bVar}:}{A list of the diagonal inner blocks (upper
      triangles only) of the positive-definite matrices on the diagonal
      of the inverse of Z'Z+W.}
    \item{\code{L}:}{The blocked sparse representation of the unit lower
      triangular matrix in the LDL' factorization of Z'Z+W.}
    \item{\code{ZZpO}:}{Object of class \code{"list"} ~~ }
   \item{\code{Omega}:}{A list of numeric matrices providing the
      components of symmetric, positive-definite matrix Omega.  Only the
      upper triangle of each component is used and stored.}
    \item{\code{D}:}{A list of the diagonal factors (upper triangle) in the
      LDL' decomposition of Z'Z+W.}
    \item{\code{REML}:}{A \code{"logical"} indicator of whether the REML
      criterion is being used to determine parameter estimates.}
    \item{\code{RXX}:}{A matrix which is the (augmented) RXX component
      or the corresponding component from the inverse of Z'Z+W}
    \item{\code{RZX}:}{A matrix which is the (augmented) RZX component
      or the corresponding component from the inverse of Z'Z+W}
    \item{\code{ZtZ}:}{The blocks in the blocked sparse symmetric
      representation of the original Z'Z matrix}
    \item{\code{XtX}:}{The original X'X matrix}
    \item{\code{ZtX}:}{The original Z'X matrix}
    \item{\code{cnames}:}{Column names of the model matrices.}
    \item{\code{devComp}:}{Numeric vector of length 4 giving the
      components used to calculate the deviance. This slot's contents
      are current if \code{status$factored} is \code{TRUE}.}
    \item{\code{deviance}:}{Numeric vector of length 2 containing the
      deviance corresponding to the maximum likelihood (ML) and REML
      criteria.  This slot's contents are current if
      \code{status$factored} is \code{TRUE}.}
    \item{\code{nc}:}{Integer vector of the number of columns in
      (augmented) model matrices and the number of observations in the
      model frame.}
    \item{\code{Gp}:}{Integer vector of 0-based row numbers of the
      \code{ZtX} and \code{RZX} matrices where new groups begin.}
    groups of rows in }
    \item{\code{status}:}{Logical vector of length 2 indicating if the
      object has been factored and if the factorization has been
      inverted.}
    \item{\code{call}:}{A copy of the matched call to the function that
      created the object.}
  }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "lmer")}: Extract the
      parameters that determine the Omega matrices.  The optional
      logical argument \code{unconstr} determines the constrained or
      unconstrained parameterization.}
    \item{coef<-}{\code{signature(object = "lmer", value = "numeric")}: Assign the
      parameters that determine the Omega matrices.  The optional
      logical argument \code{unconstr} determines the constrained or
      unconstrained parameterization. }
    \item{coerce}{\code{signature(from = "lmer", to = "dtTMatrix")}:
      Ensure that the \code{lmer} object has been factored but not
      inverted and extract the unit lower triangular matrix \eqn{L}{L}
      from the LDL decomposition of \eqn{Z^\prime
	Z+\Omega}{crossprod(Z)+Omega}.}
    \item{coerce}{\code{signature(from = "lmer", to = "dsTMatrix")}:
      Ensure that the \code{lmer} object has been factored but not
      inverted and extract the symmetric matrix \eqn{Z^\prime
	Z+\Omega}{crossprod(Z)+Omega}.}
    \item{solve}{\code{signature(a = "lmer", b = "missing")}: Invert
      the decomposed matrices.}
  }
}
\keyword{classes}

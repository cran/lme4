\name{LMEgradient}
\docType{genericFunction}
\alias{LMEgradient}
\alias{LMEgradient-methods}
\alias{LMEgradient,pdDiag,matrix,numeric-method}
\alias{LMEgradient,pdIdent,matrix,numeric-method}
\alias{LMEgradient,pdLogChol,matrix,numeric-method}
\alias{LMEgradient,pdNatural,matrix,numeric-method}
\alias{LMEgradient,lmeLevel,matrix,missing-method}
\alias{LMEgradient,lmeLevel,missing,missing-method}
\alias{LMEhessian}
\alias{LMEhessian-methods}
\alias{LMEhessian,pdDiag,matrix,array,numeric-method}
\alias{LMEhessian,pdLogChol,matrix,array,numeric-method}
\alias{LMEhessian,lmeLevel,missing,missing,missing-method}
\alias{LMEhessian,reStruct,missing,missing,missing-method}
\title{The gradient and Hessian in lme optimization-method}
\description{
  The LMEgradient and LMEhessian generic functions return the gradient
  and the Hessian of the
  log-likelihood or log-restricted-likelihood in an lme model with
  respect to the parameters of the object represented by \code{x}.
}
\usage{
LMEgradient(x, A, nlev)
LMEhessian(x, A, H, nlev)
}
\arguments{
  \item{x}{a parameterized component of an lme model, usually the
    precision matrix of an \code{lmeLevel}.  Such precision matrices
    inherit from the \code{pdMat} class.}
  \item{A}{an upper triangular matrix with the same number of columns as
    the matrix represented by \code{x}}
  \item{H}{an array of four dimensions with each dimension same as the
    number of columns of the matrix represented by \code{x}}
  \item{nlev}{integer: the number of levels of the grouping factor
    corresponding to the random-effects structure}
}
%\details{}
\value{
  LMEgradient returns
  a numeric vector of length \code{length(coef(x))}.

  LMEhessian returns
  a symmetric matrix with number of columns \code{length(coef(x))}.
}
%\references{}
\author{Douglas Bates \email{bates@stat.wisc.edu}
  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\note{}
\seealso{\code{\link{pdMat-class}}, \code{\link{lmeLevel-class}}}
\examples{
library(lme4)
data(Oxboys)
m3 <- pdLogChol(~ age)
as(m3, 'pdmatrix') <- crossprod(model.matrix(formula(m3), as(Oxboys, "data.frame")))
LMEgradient(m3, diag(2), 19)
}
\keyword{models}

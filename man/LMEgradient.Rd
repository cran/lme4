\name{LMEgradient}
\docType{genericFunction}
\alias{LMEgradient}
\alias{LMEgradient-methods}
\alias{LMEgradient,pdDiag,matrix,numeric-method}
\alias{LMEgradient,pdIdent,matrix,numeric-method}
\alias{LMEgradient,pdLogChol,matrix,numeric-method}
%\alias{LMEgradient,pdScalar,matrix,numeric-method}
\alias{LMEgradient,lmeLevel,matrix,missing-method}
\title{The gradient in lme optimization}
\description{
  This generic function returns the gradient of the log-likelihood or
  log-restricted-likelihood in an lme model with respect to the
  parameters of the object represented by \code{x}.
}
\usage{
LMEgradient(x, A, nlev)
}
\arguments{
  \item{x}{a parameterized component of an lme model, usually the
    precision matrix of an \code{lmeLevel}.  Such precision matrices
    inherit from the \code{pdMat} class.}
  \item{A}{an upper triangular matrix with the same number of columns as
    the matrix represented by \code{x}}
  \item{nlev}{integer: the number of levels of the grouping factor
    corresponding to the random-effects structure}
}
%\details{}
\value{
  a numeric vector of length \code{length(coef(x))}.
}
%\references{}
\author{Douglas Bates \email{bates@stat.wisc.edu}
  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\note{}
\seealso{\code{\link{pdMat-class}}, \code{\link{lmeLevel-class}}}
\examples{
library(lme4)
data(Oxboys, package = "nlme")
m3 <- pdLogChol(~ age)
as(m3, 'pdmatrix') <- crossprod(model.matrix(formula(m3), Oxboys))
LMEgradient(m3, diag(2), 19)
}
\keyword{models}

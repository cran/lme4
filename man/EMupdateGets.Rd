\name{EMupdateGets}
\docType{genericFunction}
\alias{EMupdate<-}
\alias{EMupdate<--methods}
\alias{EMupdate<-,pdDiag,numeric,matrix-method}
\alias{EMupdate<-,pdIdent,numeric,matrix-method}
\alias{EMupdate<-,pdLogChol,numeric,matrix-method}
%\alias{EMupdate<-,pdScalar,numeric,matrix-method}
\alias{EMupdate<-,lmeLevel,missing,matrix-method}
\title{Update an object in the EM algorithm}
\description{
  A generic replacement function used in the EM and
  ECME optimization algorithms for linear mixed-effects models.
}
\usage{
EMupdate(x, nlev) <- value
}
\arguments{
  \item{x}{an \code{lmeStruct} or \code{pdMat} object}
  \item{nlev}{integer: the number of levels of the grouping factor
    corresponding to the random-effects structure}
  \item{value}{a matrix with the triangular factor from an
    orthogonal-triangular decomposition of the modal values of the
    random-effects and their precision matrices}
}
%\details{}
\value{
  an object of the same class as \code{x}
}
%\references{}
\author{Douglas Bates \email{bates@stat.wisc.edu}
  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\note{}
%\seealso{}
\examples{
library(lme4)
data(Oxboys)
m3 <- pdLogChol(~ age)
m3
as(m3, 'pdmatrix') <- crossprod(model.matrix(formula(m3), as(Oxboys, "data.frame")))
show(m3)
EMupdate(m3, length(levels(Oxboys$Subject))) <- diag(2)
show(m3)
}
\keyword{models}

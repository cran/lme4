\name{logDet}
\docType{genericFunction}
\docType{methods}
\alias{logDet}
\alias{logDet-methods}
\alias{logDet,pdBlocked,missing-method}
\alias{logDet,pdCompSymm,missing-method}
\alias{logDet,pdDiag,missing-method}
\alias{logDet,pdIdent,missing-method}
\alias{logDet,pdLogChol,missing-method}
\alias{logDet,pdMat-method}
\alias{logDet,pdNatural,missing-method}
%\alias{logDet,pdScalar,missing-method}
\title{Logarithm of the determinant}
\description{
  Return the logarithm of the determinant of an object of a suitable
  class.
}
\usage{
logDet(object, covariate = getCovariate(object), ...)
}
\arguments{
  \item{object}{An object of a suitable class.}
  \item{covariate}{Optionally, a covariate vector}
  \item{\dots}{Some methods may have optional arguments.}
}
%\details{}
\value{
  A numeric value: the logarithm of the determinant of the matrix
  represented by \code{object}. 
}
%\references{}
%\author{Douglas Bates \email{bates@stat.wisc.edu} and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\note{}
%\seealso{}
%\examples{}
\keyword{algebra}
\keyword{array}

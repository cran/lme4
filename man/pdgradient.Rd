\name{pdgradient}
\docType{methods}
\alias{pdgradient}
\alias{pdgradient-methods}
\alias{pdgradient,pdIdent-method}
\alias{pdgradient,pdLogChol-method}
\alias{pdgradient,pdMat-method}
\title{Evaluated the pdmatrix gradient}
\description{
  This function is generic.  It returns the gradient of the positive
  definite matrix represented by \code{x} with respect to the parameters
  of \code{x}.
}
\usage{
pdgradient(x)
}
\arguments{
  \item{x}{An object of a suitable class, usually one of
  the classes that inherit from \code{pdMat}.}
}
\section{Methods}{
  \describe{
    \item{x = pdIdent}{gradient for the \code{pdIdent} class}
    \item{x = pdLogChol}{gradient for the \code{pdLogChol} class}
    \item{x = pdMat}{gradient for the \code{pdMat} class, implemented
      using \code{pdgradNumeric}}
%    \item{x = pdScalar}{gradient for the \code{pdScalar} class}
  }
}
%\details{}
\value{
  A numeric array of dimension \code{c(dim(x), length(coef(x)))} that
  is the gradient of the \code{pdmatrix} represented by \code{x} with
  respect to the parameters of \code{x}.
}
%\references{}
\author{Douglas Bates \email{bates@stat.wisc.edu}
  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\note{}
\seealso{\code{\link{pdgradNumeric}}}
\examples{
}
\keyword{array}
\keyword{algebra}
\keyword{models}
\keyword{methods}

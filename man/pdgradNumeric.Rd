\name{pdgradNumeric}
\alias{pdgradNumeric}
\title{Evaluate the pdmatrix gradient numerically}
\description{
  Evaluate the pdmatrix gradient via finite differences.  This function
  is used to validate the \code{pdgradient} methods for \code{pdMat} classes.
}
\usage{
pdgradNumeric(x)
}
\arguments{
  \item{x}{A \code{pdMat} object}
}
%\details{}
\value{
  A finite-difference approximation to the gradient of the
  \code{pdmatrix} represented by \code{x} with respect to the parameters
  of \code{x}.
}
%\references{}
\author{Douglas Bates \email{bates@stat.wisc.edu}
  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\note{}
%\seealso{}
\examples{
}
\keyword{array}
\keyword{algebra}
\keyword{models}


\name{weightedGets}
\docType{methods}
\docType{genericFunction}
\alias{weighted<-}
\alias{weighted<-,reStruct,matrix-method}
\title{Modify the weights}
\description{
  Modify the weights in an object of a suitable class. This function is
  generic.
}
\usage{
weighted(x) <- value
}
\arguments{
  \item{x}{An object of a suitable class.}
  \item{value}{weights to assign to the object \code{x}.}
}
%\details{}
\value{
  The object \code{x} with the weights updated.
}
\section{Methods}{
  \describe{
    \item{x = reStruct, value = matrix}{modify the weights in an
      \code{reStruct} object}
  }
}
%\author{Douglas Bates \email{bates@stat.wisc.edu}
%  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\note{}
\seealso{\code{\link{reStruct-class}}}
%\examples{}
\keyword{manip}
\keyword{models}

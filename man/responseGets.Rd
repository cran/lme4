\name{responseGets}
\docType{methods}
\docType{genericFunction}
\alias{response<-}
\alias{response<-,reStruct,matrix-method}
\title{Modify the response}
\description{
  Modify the reponse in an object of a suitable class. This function is
  generic.
}
\usage{
response(x) <- value
}
\arguments{
  \item{x}{An object of a suitable class.}
  \item{value}{updated response vector to assign to the object \code{x}.}
}
%\details{}
\value{
  The object \code{x} with the response updated.
}
\section{Methods}{
  \describe{
    \item{x = reStruct, value = matrix}{modify the response in an
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

\name{EMstepsGets}
\docType{genericFunction}
\alias{EMsteps<-}
\alias{EMsteps<-,reStruct,list-method}
\title{Perform EM iterations}
\description{
  The object \code{x} is updated through EM or ECME iterations.  The
  number of iterations and whether verbose output should be produced are
  specified in the list passed as \code{value}.
}
\usage{
EMsteps(x) <- value
}
\arguments{
  \item{x}{An object of a suitable class to be updated according to the
    EM or ECME algorithm.  Generally an \code{reStruct} object.}
  \item{value}{A list that contains components \code{EMiter}, a
    non-negative integer specifying the number of iterations, and
    \code{EMverbose}, a logical value that, when \code{TRUE}, requests
    verbose output.}
}
\value{
  The object \code{x}, updated.
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link{reStruct-class}}, \code{\link{lmeControl}}}
\examples{}
\keyword{models}% at least one, from doc/KEYWORDS

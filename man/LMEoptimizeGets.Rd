\name{LMEoptimizeGets}
\docType{genericFunction}
\alias{LMEoptimize<-}
\alias{LMEoptimize<-,reStruct,list-method}
\alias{LMEoptimize<-,glmm,list-method}
\title{Optimize linear mixed-effects}
\description{
  Update the parameters in an object by using an optimizer such as
  \code{optim} or \code{nlm} to maximize the log-likelihood.
}
\usage{
LMEoptimize(x) <- value
}
\arguments{
  \item{x}{An object of a suitable class, such as \code{reStruct} to be
    optimized}
  \item{value}{A list that contains components \code{msMaxIter}, a
    non-negative integer specifying the maximum number of iterations, 
    \code{msVerbose}, a logical value that, when \code{TRUE}, requests
    verbose output, \code{msTol}, a numeric tolerance for declaring
    convergence, and \code{optimizer}, a character vector whose first
    element is \code{"optim"} or \code{"nlm"}}
}
%\details{}
\value{
  The object \code{x}, updated.
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link{reStruct-class}}, \code{\link{lmeControl}}}
%\examples{}
\keyword{models}

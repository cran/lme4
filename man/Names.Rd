\name{Names}
\docType{genericFunction}
\alias{Names}
\alias{Names<-}
\alias{Names-methods}
\alias{Names<--methods}
\alias{Names,ANY-method}
\alias{Names<-,ANY-method}
\title{Extract or assign Names}
\description{
  Generic functions to extract or assign names of an object, usually an
  object of a class that inherits from the \code{pdMat} class.  This
  function is for back-compatibility only.
}
\usage{
Names(object, ...)
Names(object) <- value
}
\arguments{
  \item{object}{An object with names - usually an object that inherits
    from the \code{pdMat} class.}
  \item{value}{A character vector of the same length as \code{Names(object)}}
  \item{\dots}{Some methods for this generic function may use
    additional, optional arguments.}
}
\value{
  The names assigned to \code{object}.
}
\author{Douglas Bates \email{bates@stat.wisc.edu} and
  Saikat DebRoy \email{saikat@stat.wisc.edu}}
\seealso{\code{\link{names}}, \code{\link{names<-}}}
\examples{
}
\keyword{attribute}


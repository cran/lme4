\name{model.matrixGets}
\docType{methods}
\docType{genericFunction}
\alias{model.matrix<-}
\alias{model.matrix<-,reStruct,matrix-method}
\title{Assign the model.matrix}
\description{
  Assign the model.matrix in an object of a suitable class. This function is
  generic.
}
\usage{
model.matrix(x) <- value
}
\arguments{
  \item{x}{An object of a suitable class.}
  \item{value}{A model matrix to assign to the object \code{x}.}
}
%\details{}
\value{
  The object \code{x} with the model.matrix updated.
}
\section{Methods}{
  \describe{
    \item{x = reStruct, value = matrix}{assign the model.matrix in an
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

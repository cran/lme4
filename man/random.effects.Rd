% $Id: random.effects.Rd,v 1.1 2002/11/27 22:49:45 bates Exp $
\name{random.effects}
\docType{genericFunction}
\docType{methods}
\title{Extract Random Effects}
\alias{random.effects}
\alias{ranef}
\alias{ranef-methods}
\alias{ranef,ANY-method}
\alias{ranef,reStruct-method}
\usage{
random.effects(object, \dots)
ranef(object, \dots)
}
\description{
  A generic function to extract the random effects.
}
\arguments{
 \item{object}{an object of a class from which random effects
   estimates can be extracted.}
 \item{\dots}{some methods for this generic function require additional
   arguments.}
}
\value{
  a list of matrices.
}
%\seealso{}
\examples{
## see the method function documentation
}
\keyword{models}
\keyword{methods}

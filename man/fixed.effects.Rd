% $Id: fixed.effects.Rd,v 1.2 2002/11/29 15:52:18 bates Exp $
\name{fixed.effects}
\docType{genericFunction}
\alias{fixed.effects}
\alias{fixef}
\alias{fixef,ANY-method}
\alias{fixef,reStruct-method}
\title{Extract Fixed Effects}
\description{
  A generic function to extract the fixed effects.
}
\usage{
fixed.effects(object, \dots)
fixef(object, \dots)
}
\arguments{
 \item{object}{any fitted model object from which fixed effects
   estimates can be extracted.}
 \item{\dots}{some methods for this generic function require additional
   arguments.}
}
\value{
  will depend on the method function used; see the appropriate documentation.
}
%\seealso{}
\examples{
## see the method function documentation
}
\keyword{models}

\name{fixed.effects}
\docType{genericFunction}
\alias{fixed.effects}
\alias{fixef}
\alias{fixef,ANY-method}
\alias{fixef,ssclme-method}
\alias{fixef,lmeRep-method}
\alias{fixef,lme-method}
\title{Extract Fixed Effects}
\description{
  Extract the estimates of the fixed effects parameters from a fitted model.
}
\usage{
fixef(object, \dots)
fixed.effects(object, \dots)  # deprecated
}
\arguments{
 \item{object}{any fitted model object from which fixed effects
   estimates can be extracted.}
 \item{\dots}{optional additional arguments. Currently none are used in
   any methods.}
}
\value{
  a named numeric vector of fixed effects estimates.
}
%\seealso{}
%\examples{}
\keyword{models}

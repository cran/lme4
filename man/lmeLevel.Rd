\name{lmeLevel}
\docType{genericFunction}
\alias{lmeLevel}
\alias{lmeLevel,pdMat,factor,integer,matrix-method}
\alias{lmeLevel,missing,missing,integer,matrix-method}
\title{Construct an lmeLevel object}
\description{
  Create an \code{lmeLevel} object from its components}
}
\usage{
lmeLevel(precision, groups, columns, modelMatrix)
}
\arguments{
  \item{precision}{A \code{pdMat} object representing the relative
    precision matrix for the random effects.  This argument is missing
    in calls to generate the \code{lmeLevel} for the fixed-effects terms
    and for the response.}
  \item{groups}{The grouping factor.  This argument is missing
    in calls to generate the \code{lmeLevel} for the fixed-effects terms
    and for the response.}
  \item{columns}{integer vector: the columns in \code{modelMatrix}
    corresponding to this level of random or fixed effects.}
  \item{modelMatrix}{a model matrix: the \code{original} slot of the
    containing \code{reStruct} object}
}
\details{
  An \code{lmeLevel} object contains information on one level of
  random-effects terms or the fixed-effects terms or the response in a
  linear mixed-effects structure.
}
\value{
  An \code{lmeLevel} object.
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link{lmeLevel-class}}, \code{\link{pdMat-class}},
  \code{\link{reStruct-class}}}
%\examples{}
\keyword{models}

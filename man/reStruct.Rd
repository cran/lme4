\name{reStruct}
\docType{methods}
\docType{genericFunction}
\alias{reStruct}
\alias{reStruct,formula,list,data.frame,missing,logical,missing-method}
\alias{reStruct,formula,list,data.frame,numeric,logical,missing-method}
\title{Create an reStruct object}
\description{
  Create an \code{reStruct} object from values for its slots.
}
\usage{
reStruct(fixed, random, data, weights, REML, nextraCols)
}
\arguments{
  \item{fixed}{a formula of the form \code{response ~ terms} where
    \code{response} is an expression for the response and \code{terms}
    are one or more terms that define the fixed-effects.}
  \item{random}{a named list of one-sided formulas.  The names in the
    list evaluate to the grouping factors and the right-hand side of the
    formula is an expression for a \code{pdMat} object that itself
    contains a formula.  The order of the components of the list
    determines the nesting structure in that the second component is
    nested within the first, the third within the second, and so on.}
  \item{data}{a \code{data.frame} or \code{model.frame} in which the
    \code{fixed} and \code{random} arguments can be evaluated.}
  \item{weights}{a numeric vector of prior case weights for the
    observations. Defaults to equal weights.}
  \item{REML}{a \code{logical} value: if \code{TRUE} the REML
    criterion will be used for parameter estimation, otherwise the ML
    criterion will be used.}
  \item{nextraCols}{an \code{integer} giving the number of extra columns
    to be added to the stored matrices.}
}
%\details{}
\value{
  An \code{reStruct} object.
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link{reStruct-class}}, \code{\link{lmeLevel-class}}}
%\examples{}
\keyword{models}

\name{glmmStruct}
\docType{methods}
\docType{genericFunction}
\alias{glmmStruct}
\alias{glmmStruct,formula,list-method}
\title{Create a glmm object}
\description{
  Create a \code{glmm} object from values for its slots.
}
\usage{
glmmStruct(formula, random, family, data, nextraCols, method, ...)
}
\arguments{
  \item{formula}{a formula of the form \code{response ~ terms} where
    \code{response} is an expression for the response and \code{terms}
    are one or more terms that define the fixed-effects.}
  \item{random}{a named list of one-sided formulas.  The names in the
    list evaluate to the grouping factors and the right-hand side of the
    formula is an expression for a \code{pdMat} object that itself
    contains a formula.  The order of the components of the list
    determines the nesting structure in that the second component is
    nested within the first, the third within the second, and so on.}
  \item{family}{a description of the error distribution and link
    function to be used in the model.
    See \code{\link[base:family]{family}} for details.}
  \item{data}{a \code{data.frame} or \code{model.frame} in which the
    \code{fixed} and \code{random} arguments can be evaluated.}
  \item{nextraCols}{an \code{integer} giving the number of extra columns
    to be added to the stored matrices.}
  \item{method}{a \code{character} giving the method used to fit the
    model.}
  \item{\dots}{optional arguments passed to
    \code{\link[base:model.frame]{model.frame}}}
}
%\details{}
\value{
  A \code{glmm} object.
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link{glmm-class}}, \code{\link{lmeLevel-class}}}
%\examples{}
\keyword{models}

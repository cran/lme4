\name{formulas}
\title{Utilities for formulas}
\alias{getCovariate}
\alias{getCovariateFormula}
\alias{getGroups}
\alias{getGroups,data.frame,formula-method}
\alias{getGroups,groupedData-method}
\alias{getGroupsFormula}
\alias{getGroupsFormula,ANY-method}
\alias{getResponse}
\alias{getResponseFormula}
\alias{subFormula}
\alias{splitFormula}
\description{
  These generic functions extract sub-formulas or values of components
  from a formula/data pair.
}
\usage{
getCovariate(object, form, data)
getCovariateFormula(object)
getGroups(object, form, level, data, sep)
getGroupsFormula(object, asList, sep)
getResponse(object, form)
getResponseFormula(object)
subFormula(form, pos)
splitFormula(form, sep)
}
\arguments{
  \item{object}{a formula or an object that can provide a formula}
  \item{form}{an optional formula with a conditioning expression on its
    right hand side (i.e. an expression involving the \code{|}
    operator). Defaults to \code{formula(object)}.}
  \item{level}{a positive integer vector with the level(s) of grouping to
    be used when multiple nested levels of grouping are present. This
    argument is optional for most methods of this generic function and
    defaults to all levels of nesting.}
  \item{data}{a data frame in which to interpret the variables named in
    \code{form}. Optional for most methods.}
  \item{asList}{an optional logical value. If \code{TRUE} the returned
    value with be a list of formulas; else, if \code{FALSE} the returned
    value will be a one-sided formula. Defaults to \code{FALSE}.}
  \item{sep}{character, the separator to use between group levels when
    multiple levels are collapsed.  Defaults to \code{'/'}.}
  \item{pos}{the integer index of the component of the formula to
    return.  Defaults to \code{2}.}
}
\value{
  \code{getGroupsFormula} returns a one-sided formula, or a list of
  one-sided formulas, giving the grouping structure associated with
  \code{formula(object)}. If no conditioning expression is present in
  \code{formula(object)} a \code{NULL} value is returned.
}
%\author{Douglas Bates \email{bates@stat.wisc.edu} and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\seealso{}
\examples{
form = resp ~ cov1 + cov2 | grp1/grp2
getResponseFormula(form)
getCovariateFormula(form)
getGroupsFormula(form)
getGroupsFormula(form, asList = TRUE)
}
\keyword{models}

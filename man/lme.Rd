\name{lme}
\docType{genericFunction}
\alias{lme}
\alias{lme,ANY,ANY,formula-method}
%\alias{lme,ANY,missing,ANY-method}
\alias{lme,missing,ANY,ANY-method}
\alias{lme,formula,ANY,formula-method}
\alias{lme,formula,ANY,list-method}
\alias{lme,formula,data.frame,list-method}
\alias{lme,formula,groupedData,list-method}
\alias{lme,formula,groupedData,missing-method}
\alias{lme,missing,groupedData,ANY-method}
\alias{lme1}
\alias{lme1,ANY,ANY,formula-method}
\alias{lme1,missing,ANY,ANY-method}
\alias{lme1,formula,ANY,formula-method}
\alias{lme1,formula,ANY,list-method}
\alias{lme1,formula,data.frame,list-method}
\alias{lme1,formula,groupedData,list-method}
\alias{lme1,formula,groupedData,missing-method}
\alias{lme1,missing,groupedData,ANY-method}
\title{Fit linear mixed-effects models}
\description{
  This generic function fits a linear mixed-effects model in the
  formulation described in Laird and Ware (1982) but allowing for nested
  random effects. The within-group errors are allowed to be correlated
  and/or have unequal variances. 
}
\synopsis{
lme(formula, data, random, method = c("REML", "ML"),
    control = list(), subset, weights, na.action, offset,
    model = TRUE, x = FALSE, y = FALSE, \dots)
}
\usage{
lme(formula, data, random, \dots)
}
\arguments{
  \item{formula}{a two-sided linear formula object describing the
    fixed-effects part of the model, with the response on the left of a
    \code{~} operator and the terms, separated by \code{+} operators, on
    the right.}
  \item{data}{an optional data frame containing the variables named in
    \code{formula}, \code{random}, \code{correlation}, \code{weights}, and
    \code{subset}.  By default the variables are taken from the
    environment from which \code{lme} is called.}
 \item{random}{optionally, any of the following: (i) a one-sided formula
   of the form \code{~x1+...+xn | g1/.../gm}, with \code{x1+...+xn}
   specifying the model for the random effects and \code{g1/.../gm} the
   grouping structure (\code{m} may be equal to 1, in which case no
   \code{/} is required). The random effects formula will be repeated
   for all levels of grouping, in the case of multiple levels of
   grouping; (ii) a list of one-sided formulas of the form
   \code{~x1+...+xn | g}, with possibly different random effects models
   for each grouping factor; or (iii) a named list of formulas.}
 \item{\dots}{Optional arguments for methods.  Currently none are used.}
}
\details{
  A new version of the function called \code{lme1} is in the development
  versions of this package.  The \code{lme1} function has the same
  argument sequence as \code{lme} but produces an object with an
  \code{\link[Matrix]{lmeRep}} representation of the model.  The object
  produced by \code{lme} uses the \code{\link[Matrix]{ssclme}}
  representation.
  
  Additional standard arguments to model-fitting functions can be passed
  to \code{lme} or \code{lme1}.
  \describe{
    \item{correlation}{an optional \code{corStruct} object describing
      the within-group correlation structure. See the documentation of
      \code{corClasses} for a description of the available \code{corStruct}
      classes. Defaults to \code{NULL}, corresponding to no within-group
      correlations.} 
    \item{weights}{an optional \code{varFunc} object or one-sided formula
      describing the within-group heteroscedasticity structure. If given as
      a formula, it is used as the argument to \code{varFixed},
      corresponding to fixed variance weights. See the documentation on
      \code{varClasses} for a description of the available \code{varFunc}
      classes. Defaults to \code{NULL}, corresponding to homocesdatic
      within-group errors.} 
    \item{subset}{an optional expression indicating the subset of the rows of
      \code{data} that should be used in the fit. This can be a logical
      vector, or a numeric vector indicating which observation numbers are
      to be included, or a  character  vector of the row names to be
      included.  All observations are included by default.}
    \item{method}{a character string.  If \code{"REML"} the model is fit by
      maximizing the restricted log-likelihood.  If \code{"ML"} the
      log-likelihood is maximized.  Defaults to \code{"REML"}.}
    \item{na.action}{a function that indicates what should happen when the
      data contain \code{NA}s.  The default action (\code{na.fail}) causes
      \code{lme} to print an error message and terminate if there are any
      incomplete observations.}
    \item{control}{a list of control values for the estimation algorithm to
      replace the default values returned by the function \code{lmeControl}.
      Defaults to an empty list.}
    \item{model, x}{logicals.  If \code{TRUE} the corresponding
      components of the fit (the model frame, the model matrices)
      are returned.}
  }
}
\value{
  An \code{\link{lme-class}{lme}} object.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{}
%\note{}
\seealso{\code{\link{lme-class}}, \code{\link[stats]{lm}}}
\examples{
data(bdf)
fm <- lme(langPOST ~ IQ.ver.cen + avg.IQ.ver.cen, data = bdf,
          random = ~ IQ.ver.cen | schoolNR)
summary(fm)
}
\keyword{models}
\keyword{methods}

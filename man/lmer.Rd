\name{lmer}
\docType{genericFunction}
\alias{lmer}
\alias{lmer,formula-method}
\title{Fit linear mixed-effects models}
\description{
  This generic function fits a linear mixed-effects model with nested or
  crossed grouping factors for the random effects.
}
\synopsis{
lmer(formula, data, method = c("REML", "ML"),
    control = list(), subset, weights, na.action, offset,
    model = TRUE, x = FALSE, y = FALSE, \dots)
}
\usage{
lmer(formula, data, \dots)
}
\arguments{
  \item{formula}{a two-sided linear formula object describing the
    fixed-effects part of the model, with the response on the left of a
    \code{~} operator and the terms, separated by \code{+} operators, on
    the right.  The vertical bar character \code{"|"} separates an
    expression for a model matrix and a grouping factor.}
  \item{data}{an optional data frame containing the variables named in
    \code{formula}.  By default the variables are taken from the
    environment from which \code{lmer} is called.}
 \item{\dots}{Optional arguments for methods.  Currently none are used.}
}
\details{
  This is a revised version of the \code{lme} function from the
  \code{nlme} package.
  
  Additional standard arguments to model-fitting functions can be passed
  to \code{lmer}.
  \describe{
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
      replace the default values returned by the function \code{lmerControl}.
      Defaults to an empty list.}
    \item{model, x}{logicals.  If \code{TRUE} the corresponding
      components of the fit (the model frame, the model matrices)
      are returned.}
  }
}
\value{
  An \code{\link[Matrix]{lme-class}{lmer}} object.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{}
%\note{}
\seealso{\code{\link[Matrix]{lmer-class}}, \code{\link[stats]{lm}}}
\examples{
(fm <- lmer(langPOST ~
           IQ.ver.cen + avg.IQ.ver.cen + (IQ.ver.cen | schoolNR),
           data = bdf))
}
\keyword{models}
\keyword{methods}

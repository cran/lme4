\name{lmer}
\docType{genericFunction}
\alias{lmer}
\alias{lmer,formula,ANY,ANY-method}
\alias{lmer,formula,ANY,missing-method}
\title{Fit linear mixed-effects models}
\description{
  This generic function fits a linear mixed-effects model with nested or
  crossed grouping factors for the random effects.
}
\synopsis{
lmer(formula, data, family, method = c("REML", "ML", "PQL", "Laplace", "AGQ"),
    control = list(), subset, weights, na.action, offset,
    model = TRUE, x = FALSE, y = FALSE, \dots)
}
\usage{
lmer(formula, data, family,
     method = c("REML", "ML", "PQL", "Laplace", "AGQ"),\dots)
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
  \item{family}{
    a GLM family, see \code{\link[stats]{glm}}.  If \code{family} is
    missing then a linear mixed model is fit; otherwise a generalized
    linear mixed model is fit.}
  \item{method}{a character string.  For a linear mixed model the
    default is \code{"REML"} indicating that the model should be fit by
    maximizing the restricted log-likelihood.  The alternative is
    \code{"ML"} indicating that the log-likelihood should be
    maximized.  (This method is sometimes called "full" maximum
    likelihood.)  For a generalized linear mixed model the criterion
    is always the log-likelihood but this criterion does not have a
    closed form expression and must be approximated.  The default approximation is
    \code{"PQL"} or penalized quasi-likelihood.  Alternatives are
    \code{"Laplace"} or \code{"AGQ"} indicating the Laplacian and
    adaptive Gaussian quadrature approximations respectively.  The
    \code{"PQL"} method is fastest but least accurate.  The
    \code{"Laplace"} method is intermediate in speed and accuracy.
    The \code{"AGQ"} method is the most accurate but can be
    considerably slower than the others.}
  \item{\dots}{Optional arguments for methods.  Currently none are used.}
  }
}
\details{
  This is a revised version of the \code{lme} function from the
  \code{nlme} package.  This version uses a different method of
  specifying random-effects terms and allows for fitting generalized
  linear mixed models as well as linear mixed models.
  
  Additional standard arguments to model-fitting functions can be passed
  to \code{lmer}.
  \describe{
    \item{subset}{an optional expression indicating the subset of the rows of
      \code{data} that should be used in the fit. This can be a logical
      vector, or a numeric vector indicating which observation numbers are
      to be included, or a  character  vector of the row names to be
      included.  All observations are included by default.}
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
(fm1 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
             OrchardSprays))
}
\keyword{models}
\keyword{methods}

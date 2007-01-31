\name{VarCorr}
\docType{genericFunction}
\alias{VarCorr}
\alias{VarCorr,mer-method}
\alias{VarCorr,mer2-method}
\title{Extract variance and correlation components}
\usage{
\S4method{VarCorr}{mer}(x, REML = x@status["REML"], \dots)
\S4method{VarCorr}{mer2}(x, REML = NULL, \dots)
}
\arguments{
  \item{x}{a fitted model object, usually an object inheriting from
    class \code{\linkS4class{lmer}} (or, temporarily, class
    \code{\linkS4class{lmer2}}).
  }
  \item{REML}{logical indicating if REML should be used.}
  \item{\dots}{Additional, optional arguments for some methods.  At
    present none are used.}
}
\description{
  Extract the estimated variances, standard
  deviations, and correlations of the random-effects terms in a
  linear mixed-effects model, of class \code{lme}, or a generalized
  linear mixed-effects model.  When appropriate,
  the within-group error variance and standard deviation are also
  calculated.
}
\value{
  an object of class \code{VarCorr}.
}
\seealso{
  the \code{\link{lmer}} function and \code{\linkS4class{lmer}} class;
  the result class \code{\linkS4class{VarCorr}}.
}
\examples{
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
             data = sleepstudy))
(VC <- VarCorr(fm2))
}
\keyword{models}

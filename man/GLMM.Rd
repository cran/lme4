\name{GLMM}
\alias{GLMM}
\alias{GLMM,ANY,ANY,missing,ANY-method}
\alias{GLMM,missing,ANY,groupedData,ANY-method}
\alias{GLMM,formula,ANY,groupedData,missing-method}
\alias{GLMM,formula,ANY,ANY,formula-method}
\alias{GLMM,formula,ANY,ANY,list-method}
\title{
Fit Generalized Linear Mixed Models via PQL
}
\description{
Fit a GLMM model with multivariate normal random effects, using
Penalized Quasi-Likelihood.
}
\usage{
GLMM(formula, family, data, random, control, niter, method, verbose, \dots)
}
\arguments{
  \item{formula}{
    a two-sided linear model formula giving fixed-effects part of the model.
  }
  \item{family}{
    a GLM family, see \code{\link{glm}}.
  }
  \item{random}{
    A formula or list of formulae describing the random effects.
  }
  \item{data}{
    an optional data frame used as the first place to find
    variables in the formulae.
  }
  \item{control}{
    an optional argument to be passed to \code{lme}.
  }
  \item{niter}{
    maximum number of PQL iterations.  Default is 20.
  }
  \item{method}{
    character: Estimation method to be used.  Possible values are
    \code{"PQL"}, the default, or \code{"Laplace"}.
    \code{"PQL"} provides penalized quasi-likelihood estimates.
    \code{"Laplace"} provides PQL
    estimation followed by optimization of the second-order Laplacian
    approximation to the marginal log-likelihood.
  }
  \item{verbose}{
    logical: print out record of iterations?  Default is \code{FALSE}.
  }
  \item{\dots}{
    Optional further arguments such as \code{subset} and \code{na.action}.
  }
}
\value{
A object of class \code{"glmm"}: see \code{\link{reStruct-class}}.
}
\references{
Schall, R. (1991) Estimation in generalized linear models with
random effects.
\emph{Biometrika}
\bold{78}, 719--727.


Breslow, N. E. and Clayton, D. G. (1993) Approximate inference in
generalized linear mixed models.
\emph{Journal of the American Statistical Association}
\bold{88}, 9--25.


Wolfinger, R. and O'Connell, M. (1993) Generalized linear mixed models: a
pseudo-likelihood approach.
\emph{Journal of Statistical Computation and Simulation}
\bold{48}, 233--243.
}
\seealso{
  \code{\link{lme}}
}
\examples{
library(lme4)
data(guImmun)
fm1 = GLMM(immun ~ kid2p + mom25p + ord + ethn +
                  momEd + husEd + momWork + rural + pcInd81,
          data = guImmun, family = binomial,
          random = ~1|comm/mom)
summary(fm1)
fm2 = GLMM(immun ~ kid2p + mom25p + ord + ethn +
                  momEd + husEd + momWork + rural + pcInd81,
          data = guImmun, family = binomial,
          random = ~1|comm, method = 'Laplace')
summary(fm2)
}
\keyword{models}

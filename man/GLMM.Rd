\name{GLMM}
\alias{GLMM}
\alias{GLMM,ANY,ANY,ANY,formula-method}
\alias{GLMM,formula,ANY,ANY,list-method}
\alias{GLMM,formula,ANY,groupedData,list-method}
\alias{GLMM,formula,ANY,groupedData,missing-method}
\alias{GLMM,missing,ANY,ANY,ANY-method}
\title{
Fit Generalized Linear Mixed Models via PQL
}
\description{
Fit a GLMM model with multivariate normal random effects, using
Penalized Quasi-Likelihood.
}
\synopsis{
GLMM(formula, family, data, random, method = c("PQL", "Laplace"),
     control = list(), subset, weights, na.action, offset,
     model = TRUE, x = FALSE, y = FALSE, \dots)
}
\usage{
GLMM(formula, family, data, random, \dots)
}
\arguments{
  \item{formula}{
    a two-sided linear model formula giving fixed-effects part of the model.
  }
  \item{family}{
    a GLM family, see \code{\link[stats]{glm}}.
  }
  \item{data}{
    an optional data frame used as the first place to find
    variables in the formulae.
  }
  \item{random}{
    A formula or named list of formulae describing the random effects.
  }
  \item{\dots}{
    Optional further arguments such as \code{subset} and \code{na.action}.
  }
}
\details{
  Additional arguments, some of them standard in model-fitting
  functions, can be passed to \code{GLMM}.
  \describe{
    \item{subset}{an optional expression indicating the subset of the rows of
      \code{data} that should be used in the fit. This can be a logical
      vector, or a numeric vector indicating which observation numbers are
      to be included, or a  character  vector of the row names to be
      included.  All observations are included by default.}
    %  \item{niter}{
    %    maximum number of PQL iterations.  Default is 20.
    %  }
    \item{method}{
      character: Estimation method to be used.  Possible values are
      \code{"PQL"}, the default, or \code{"Laplace"}.
      \code{"PQL"} provides penalized quasi-likelihood estimates.
      \code{"Laplace"} provides PQL
      estimation followed by optimization of the second-order Laplacian
      approximation to the marginal log-likelihood.
    }
    \item{na.action}{a function that indicates what should happen when the
      data contain \code{NA}s.  The default action (\code{na.fail}) causes
      \code{lme} to print an error message and terminate if there are any
      incomplete observations.}
    \item{control}{
      an optional argument to be passed to \code{lme}.
    }
    \item{model, x}{logicals.  If \code{TRUE} the corresponding
      components of the fit (the model frame, the model matrices)
      are returned.}
  }
}
\value{
  An object of class \code{"lme"}: see \code{\link[Matrix]{ssclme-class}}.
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
(fm1 <- GLMM(decrease ~ treatment, family = poisson(),
             data = OrchardSprays,
             random = list(rowpos = ~ 1, colpos = ~1)))
}
\keyword{models}

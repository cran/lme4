\name{getME}
\alias{getME}
\alias{isREML}
\alias{isLMM}
\alias{isGLMM}
\alias{isNLMM}
\title{Extract or Get Generalize Components from a Fitted Mixed Effects Model}
\description{
  Extract (or \dQuote{get}) \dQuote{components} -- in a generalized
  sense -- from a fitted mixed    effects model, i.e. (in this version
  of the package) from an \R object of class \code{"\linkS4class{mer}"}.

  The goal is to provide \dQuote{everything a user may want} from a
  fitted \code{"mer"} object \emph{as far} as it is not available by
  methods, such as \code{\link{fixef}}, \code{\link{ranef}},
  \code{\link{vcov}}, etc.
}
\usage{
getME(object, name)

isREML(x, ...)

isLMM(x, ...)
isNLMM(x, ...)
isGLMM(x, ...)
}
\arguments{
  \item{object,x}{a fitted mixed-effects model of class
    \code{"\linkS4class{mer}"}, i.e. typically the result of
    \code{\link{lmer}()}, \code{\link{glmer}()} or \code{\link{nlmer}()}.
  }
  \item{name}{a character string specifying the name of the
    \dQuote{component}.  Note this may not be the name of
    \code{\link{slot}} of \code{object}. }
  \item{...}{additional, optional arguments.  (None are used in the
  \code{mer} methods.)}
}
\value{
  Unspecified, as very much depending on the \code{\link{name}}.
}
\seealso{
  \code{\link{getCall}()} (in \R >= 2.14.0; otherwise in \pkg{MatrixModels}).

  More standard methods for \code{*mer()} objects, such as
  \code{\link{ranef}},
  \code{\link{fixef}},
  \code{\link{vcov}}, etc.
}
\note{\code{isREML(.)} is a synonym for \code{getME(.,"is_REML")}}
\examples{
## shows many methds you should consider *before* getME():
showMethods(class = "mer")

(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
              Orange, start = c(Asym = 200, xmid = 725, scal = 350))

Z <- getME(fm1, "Z")
stopifnot(is(Z, "CsparseMatrix"),
          c(180,36) == dim(Z),
	  all.equal(fixef(fm1), getME(fm1, "beta"),
		    check.attr=FALSE, tol = 0))

## All that can be accessed [potentially ..]:
(nmME <- eval(formals(getME)$name))

\dontshow{
## internal consistency check ensuring that all work:
## "try(.)" because some are not yet implemented:
str(fm.parts <- sapply(nmME, function(nm) try(getME(fm1, nm)),
                    simplify=FALSE))

str(gm.parts <- sapply(nmME, function(nm) try(getME(gm1, nm)),
                    simplify=FALSE))

str(nm.parts <- sapply(nmME, function(nm) try(getME(nm1, nm)),
                    simplify=FALSE))
}% dont..

isLMM(fm1)
isGLMM(gm1)
## check all :
is.MM <- function(x) c(LMM = isLMM(x), GLMM= isGLMM(x), NLMM= isNLMM(x))
stopifnot(cbind(is.MM(fm1), is.MM(gm1), is.MM(nm1))
	  == diag(rep(TRUE,3)))
}
\keyword{utilities}

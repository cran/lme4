\name{nlme-internal}
%\alias{allVarsRec}
%\alias{createConLin}
\alias{getFixDF}
\alias{glmmLa2LogLikComp}
\alias{glmmLa2RespWt}
\alias{glmmLa2Wt2}
%\alias{getParsGnls}
%\alias{getParsNlme}
%\alias{glsApVar}
%\alias{glsEstimate}
%\alias{gnlsApVar}
%\alias{lmeApVar}
%\alias{nonlinModel}
%\alias{print.correlation}
\title{Internal nlme functions}
\description{
  Internal nlme functions.
}
\usage{
%allVarsRec(object)
%createConLin(fixed,data,random,\dots)
getFixDF(object)
glmmLa2LogLikComp(x, eta)
glmmLa2RespWt(fam, eta, origy, w, off)
glmmLa2Wt2(fam, eta, w)
%getParsGnls(plist, pmap, beta, N)
%getParsNlme(plist, fmap, rmapRel, bmap, groups, beta, bvec, b, level, N)
%glsApVar(glsSt, sigma, conLin, .relStep, minAbsPar, natural)
%glsEstimate(object, conLin, control)
%gnlsApVar(gnlsSt, lsigma, conLin, .relStep, minAbsPar, natural)
%lmeApVar(lmeSt, sigma, conLin, .relStep, minAbsPar, natural)
%nonlinModel(modelExpression, env, paramNames)
%\method{print}{correlation}(x, title, rdig, \dots)
}
\arguments{
  \item{object}{An object of a suitable class.}
  \item{x}{An object of class glmm}
  \item{fam}{A family object}
  \item{eta}{The linear predictors}
  \\item{origy}{The original response}
  \item{w}{The initial weights}
  \item{off}{The initial offset}
}
\details{
  These functions are internal to the package and are not intended to be
  called by the user. They will be hidden by a namespace when namespaces
  and S4 classes can work together.}
\keyword{models}

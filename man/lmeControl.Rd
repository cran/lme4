% $Id: lmeControl.Rd,v 1.5 2003/09/20 11:54:43 saikat Exp $
\name{lmeControl}
\alias{lmeControl}
\title{Control values for lme}
\description{
  The values supplied in the function call replace the defaults and a
  list with all possible arguments is returned. The returned list is
  used as the \code{control} argument in the \code{lme} function.
}
\synopsis{
lmeControl(maxIter = 50,
           msMaxIter = 50,
           tolerance = sqrt((.Machine$double.eps)),
           niterEM = 20,
           msTol = sqrt(.Machine$double.eps),
           msVerbose = getOption("verbose"),
           PQLmaxIt = 20,
           .relStep = (.Machine$double.eps)^(1/3), 
           nlmStepMax = NULL,
           EMverbose = getOption("verbose"),
           analyticGradient = TRUE,
           analyticHessian=FALSE))
}
\usage{
lmeControl(maxIter, msMaxIter, tolerance, niterEM, msTol,
           msVerbose, PQLmaxIt, .relStep, nlmStepMax, 
           EMverbose, analyticGradient, analyticHessian)
}
\arguments{
 \item{maxIter}{maximum number of iterations for the \code{lme}
   optimization algorithm. Default is 50.} 
 \item{msMaxIter}{maximum number of iterations
   for the \code{nlm} optimization step inside the \code{lme}
   optimization. Default is 50.}
 \item{tolerance}{tolerance for the convergence criterion in the
   \code{lme} algorithm. Default is 1e-6.}
 \item{niterEM}{number of iterations for the EM algorithm used to refine
   the initial estimates of the random effects variance-covariance
   coefficients. Default is 25.}
 \item{msTol}{tolerance for the convergence criterion in \code{nlm},
   passed as the \code{rel.tolerance} argument to the function (see
   documentation on \code{nlm}). Default is 1e-7. }
 \item{msVerbose}{a logical value passed as the \code{trace} argument to
   \code{nlm} (see documentation on that function). Default is
   \code{getOption("verbose")}.} 
 \item{PQLmaxIt}{maximum number of iterations for the PQL algorithm in
   \code{GLMM}.}
 \item{.relStep}{relative step for numerical derivatives
   calculations. Default is \code{.Machine$double.eps^(1/3)}.}
 \item{nlmStepMax}{stepmax value to be passed to nlm. See
     \code{\link{nlm}} for details. Default is 100.0}
 \item{EMverbose}{a logical value indicating if verbose output should be
   produced during the EM iterations.  Default is \code{getOption("verbose"}.}
 \item{analyticGradient}{a logical value indicating if the analytic
   gradient of the objective should be used.  This option is for testing
   purposes and would not normally be changed from the default.  Default
   is \code{TRUE}.}
 \item{analyticHessian}{a logical value indicating if the analytic
   hessian of the objective should be calculated. This is an
   experimental feature and at present the default is \code{FALSE}. In
   future we may use the analytic Hessian in the optimization.}
}
\value{
  a list with a component for each of the possible arguments.
}
\seealso{\code{\link{lme}}, \code{\link{nlm}}, \code{\link{optim}},
  \code{\link{lmeScale}}, \code{\link{GLMM}}}
\examples{
# decrease the maximum number iterations in the ms call and
# request that information on the evolution of the ms iterations be printed
str(lmeControl(msMaxIter = 20, msVerbose = TRUE))
}
\keyword{models}

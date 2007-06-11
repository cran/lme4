\name{nlmer}
\alias{nlmer}
\title{Fit Nonlinear Mixed Models}
\description{
  Fit a nonlinear mixed model with nested or
  crossed grouping factors for the random effects.
}
\usage{
nlmer(formula, data, control, start, verbose,
     subset, weights, na.action, contrasts,
     model, \dots)
}
\arguments{
  \item{formula}{a three part formula object describing the response,
    the nonlinear model and the fixed and random effects in the model.
    In the fixed and random effects specification the vertical bar
    character \code{"|"} separates an expression for a model matrix and
    a grouping factor.  At present evaluation of the nonlinear model
    function must return a gradient attribute.}
  \item{data}{an optional data frame containing the variables named in
    \code{formula}.  By default the variables are taken from the
    environment from which \code{nlmer} is called.}
  \item{control}{a list of control parameters. See below for details.}
  \item{start}{a named numeric vector of starting values for the fixed
    effects parameters or a list that contains an element called
    \code{"fixef"} that has this form.  Optionally the list may contain
    an element called \code{"ST"} providing a starting value for the
    \code{ST} slot.}
  \item{verbose}{logical scalar - \code{TRUE} indicates verbose output
    from the iterations during the optimization process (highly
    recommended when difficulties are encountered).  Default is \code{FALSE}.}
  \item{subset, weights, na.action, contrasts}{further model
    specification arguments as in \code{\link[stats]{lm}}; see there for
    details.}
  \item{model}{logical indicating if the model component
    should be returned (in slot \code{frame}).}
  \item{\dots}{potentially further arguments for methods.  Currently
    none are used.}
}
\details{
  This is a revised version of the \code{\link[nlme]{nlme}} function from the
  \pkg{nlme} package.  This version uses a different method of
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
    \item{control}{a named list of control parameters for the estimation
      algorithm, specifying only the ones to be changed from their
      default values.  Hence defaults to an empty list.\cr
      Possible control options and their default values are:
      \describe{
	\item{\code{maxIter}:}{maximum number of iterations for the \code{lme}
	  optimization algorithm.  Default is 50.}
	\item{\code{tolerance}:}{relative convergence tolerance for the
	  \code{lme} optimization algorithm.  Default is
	  \code{sqrt(.Machine$double.eps)}.}
	\item{\code{msMaxIter}:}{maximum number of iterations
	  for the \code{nlminb} optimization step inside the \code{lme}
	  optimization.  Default is 200.}
	\item{\code{msVerbose}:}{a logical value passed as the
	  \code{trace} argument to \code{nlminb} (see documentation on
	  that function).  Default is \code{getOption("verbose")}.}
	\item{\code{niterEM}:}{number of iterations for the EM algorithm
	  used to refine the initial estimates of the random effects
	  variance-covariance coefficients.  Default is 15.}
	\item{\code{EMverbose}:}{a logical value indicating if verbose output
	  should be produced during the EM iterations.  Default is
	  \code{getOption("verbose")}.}

	\item{\code{PQLmaxIt}:}{maximum number of iterations for the PQL
	  algorithm when fitting generalized linear mixed models.
	  Default is 30.}%% unused currently?
	\item{\code{usePQL}:}{Should the PQL method be used before
	  switching to general optimization when fitting generalized linear mixed
	  models using \code{method = "Laplace"}?  Default is \code{FALSE}.}
      }
    }
    \item{model, x}{logicals.  If \code{TRUE} the corresponding
      components of the fit (the model frame, the model matrices)
      are returned.}
  }
}
\value{
  An object of class \code{"\linkS4class{nlmer}"}.
  There are many methods applicable to \code{"nlmer"} objects, see the
  above help page. %% Hmm.. or \code{\link{showMethods}(class = "lmer")}
}
%\references{ ~put references to the literature/web site here ~ }
%\author{}
%\note{}
\seealso{
  The \code{\linkS4class{nlmer}} class
}
\examples{
(fm1 <- lme4:::nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
                    Orange, verb = 1,
                    start = c(Asym = 192.6872, xmid = 728.7544, scal = 353.5320)))
(fm2 <- lme4:::nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~ (lKe+lKa+lCl|Subject),
                    Theoph, start = c(lKe = -2.5, lKa = 0.5, lCl = -3), verb = 1))
(fm3 <- lme4:::nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
                     (lKe|Subject) + (lKa|Subject) + (lCl|Subject), Theoph,
                     start = c(lKe = -2.5, lKa = 0.5, lCl = -3), verb = 1))
(fm4 <- lme4:::nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
                     (lKa+lCl|Subject), Theoph,
                     start = c(lKe = -2.5, lKa = 0.5, lCl = -3), verb = 1))
(fm5 <- lme4:::nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
                    (lKa|Subject) + (lCl|Subject), Theoph,
                    start = c(lKe = -2.5, lKa = 0.5, lCl = -3), verb = 1))
}
\keyword{models}
\keyword{methods}

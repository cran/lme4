\name{VarCorr}
\docType{genericFunction}
\alias{VarCorr}
\title{Extract variance and correlation components}
% \arguments{
%  \item{x}{a fitted model object, usually an object inheriting from
%      class \code{lme}. 
%  }
% }
\description{
  This function calculates the estimated variances, standard
  deviations, and correlations between the random-effects terms in a
  linear mixed-effects model, of class \code{lme}, or a generalized
  linear mixed-effects model, of class \code{glmm}. When appropriate,
  the within-group error variance and standard deviation are also
  calculated.
}
\value{
  an object of class \code{VarCorr}.
}
\author{Saikat DebRoy \email{saikat@stat.wisc.edu} and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{lme}}, \code{\link{GLMM}}, \code{\link{VarCorr-class}}}
\examples{
data(bdf)
fm <- lme(langPOST ~ IQ.ver.cen + avg.IQ.ver.cen, data = bdf,
          random = ~ IQ.ver.cen | schoolNR)
VarCorr(fm)
}
\keyword{models}

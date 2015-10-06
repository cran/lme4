\name{plot.lmList4}
\alias{plot.lmList4}
\alias{plot.lmList4.confint}
\title{plots for lmList4 objects}
\usage{
\method{plot}{lmList4}(x, form, abline, id, idLabels, grid, \dots)
\method{plot}{lmList4.confint}(x, y, order, ...)
}
\arguments{
 \item{x}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.}
 \item{y}{ignored: for agreement with generic method}
 \item{form}{an optional formula specifying the desired type of
   plot. Any variable present in the original data frame used to obtain
   \code{x} can be referenced. In addition, \code{x} itself
   can be referenced in the formula using the symbol
   \code{"."}. Conditional expressions on the right of a \code{|}
   operator can be used to define separate panels in a Trellis
   display. Default is \code{resid(., type = "pool") ~ fitted(.) },
   corresponding to a plot of the standardized residuals (using a pooled
   estimate for the residual standard error) versus fitted values.
 } 
 \item{abline}{an optional numeric value, or numeric vector of length
   two. If given as a single value, a horizontal line will be added to the
   plot at that coordinate; else, if given as a vector, its values are
   used as the intercept and slope for a line added to the plot. If
   missing, no lines are added to the plot.
 }
 \item{id}{an optional numeric value, or one-sided formula. If given as
   a value, it is used as a significance level for a two-sided outlier
   test for the standardized residuals. Observations with
   absolute standardized residuals greater than the \eqn{1 - value/2}
   quantile of the standard normal distribution are identified in the
   plot using \code{idLabels}. If given as a one-sided formula, its
   right hand side must evaluate to a  logical, integer, or character
   vector which is used to identify observations in the plot. If
   missing, no observations are identified.
 }
 \item{idLabels}{an optional vector, or one-sided formula. If given as a
   vector, it is converted to character and used to label the
   observations identified according to \code{id}. If given as a
   one-sided formula, its right hand side must evaluate to a vector
   which is converted to character and used to label the identified
   observations. Default is \code{getGroups(x)}. 
 }
 \item{grid}{an optional logical value indicating whether a grid should
   be added to plot. Default depends on the type of Trellis plot used:
   if \code{xyplot} defaults to \code{TRUE}, else defaults to
   \code{FALSE}.
 }
 \item{order}{which coefficient to order the results by}
 \item{\dots}{optional arguments passed to the Trellis plot function.}
}
\description{
  diagnostic and confidence-interval plots for \code{lmList} fits
}
\details{
  The \code{plot} method for \code{lmList4} objects is
  copied from that for \code{lmList} objects; see
  \code{\link[nlme]{plot.lmList}} for details.

}
\examples{
fm.plm  <- lmList(Reaction ~ Days | Subject, sleepstudy)
## diagnostic plot: standardized residuals vs. fitted
plot(fm.plm,id=0.05)
ci <- confint(fm.plm)
## plot CIs, ordered by slope (coefficient 2)
plot(ci,order=2,ylab="Subject")
}
\author{
  Original versions in \code{nlme} package by Jose Pinheiro
  and Douglas Bates
}


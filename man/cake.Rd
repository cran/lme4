\name{cake}
\alias{cake}
\docType{data}
\title{Breakage angle of chocolate cakes}
\description{
  Data on the breakage angle of chocolate cakes made with different
  recipes and baked at different temperatures.
}
\usage{data(cake)}
\format{
  A data frame with 270 observations on the following 5 variables.
  \describe{
    \item{\code{replicate}}{a factor with levels \code{1} to \code{15}}
    \item{\code{batch}}{a factor with levels \code{1}, \code{2} and \code{3}}
    \item{\code{recipe}}{a factor with levels \code{1}, \code{2} and \code{3}}
    \item{\code{temperature}}{an ordered factor with levels \code{175} < \code{185} < \code{195} < \code{205} < \code{215} < \code{225}}
    \item{\code{angle}}{a numeric vector giving the angle at which the
      cake broke.}
  }
}
\details{
  The \code{replicate} factor is nested within the \code{branch} factor.
}
\source{
  Original data given in Cochran and Cox (1957).
}
\references{
Also cited in Lee, Nelder and Pawitan (2006)
}
\examples{
str(cake)
print(fm1 <- lmer(angle ~ recipe * temperature + (1|replicate/batch), cake,
             method = "ML"), corr = FALSE)
print(fm2 <- lmer(angle ~ recipe + temperature + (1|replicate/batch), cake,
             method = "ML"), corr = FALSE)
print(fm3 <- lmer(angle ~ recipe + as.numeric(temperature) +
             (1|replicate/batch), cake, method = "ML"), corr = FALSE)
anova(fm3, fm2, fm1)
}
\keyword{datasets}

\name{Mmmec}
\alias{Mmmec}
\docType{data}
\title{Malignant melanoma deaths in Europe}
\description{
  Malignant Melanoma Mortality in the European Community associated with
  the impact of UV radiation exposure.
}
\usage{data(Mmmec)}
\format{
  A data frame with 354 observations on the following 6 variables.
  \describe{
    \item{nation}{a factor with levels \code{Belgium}, \code{W.Germany},
      \code{Denmark}, \code{France}, \code{UK}, \code{Italy}, \code{Ireland},
      \code{Luxembourg}, and \code{Netherlands}}
    \item{region}{Region ID - a factor.}
    \item{county}{County ID - a factor.}
    \item{deaths}{Number of male deaths due to MM during 1971--1980}
    \item{expected}{Number of expected deaths.}
    \item{uvb}{Centered measure of the UVB dose reaching the earth's
      surface in each county.}
  }
}
%\details{}
\source{
  \url{http://multilevel.ioe.ac.uk/softrev/mmmec.html}
}
\references{
  Langford, I.H., Bentham, G. and McDonald, A. 1998: Multilevel
  modelling of geographically aggregated health data: a case study on
  malignant melanoma mortality and UV exposure in the European
  community. \emph{Statistics in Medicine} 17: 41-58.
}
\examples{
data(Mmmec)
str(Mmmec)
summary(Mmmec)
fm1 <- GLMM(deaths ~ offset(log(expected)) + uvb,
 family = poisson, data = Mmmec, random = ~ 1 | region,
 niter = 25)
summary(fm1)
}
\keyword{datasets}

\name{groupedData}
\alias{groupedData}
\alias{convertFromS3groupedData}
\title{Creation of groupedData objects}
\description{
  \code{groupedData} creates groupedData objects.
  \code{convertFromS3groupedData} converts older S3 style objects of the
  same class name used in the \code{nlme} package.
}
\usage{
groupedData(formula, data,
            order.groups = TRUE,
            FUN = function(x, ...) max(x, na.rm = TRUE),
            outer = NULL, inner = NULL,
            labels = list(), units = list(), ...)
convertFromS3groupedData(from)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{}
  \item{data}{}
  \item{order.groups}{}
  \item{FUN}{}
  \item{outer}{}
  \item{inner}{}
  \item{labels}{}
  \item{units}{}
  \item{\dots}{}
  \item{from}{an object of class "groupedData" as used in package \code{nlme}}
}
\details{
  See \code{nlme} docs for now. \code{order.groups, FUN} are currently
  ignored.
}
\value{
  An object of formal class \code{groupedData}
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }

\seealso{\code{\link{groupedData-class}}}
\examples{}
\keyword{manip}


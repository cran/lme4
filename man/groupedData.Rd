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
groupedData(formula, data = sys.parent(1),
            order.groups = TRUE,
            FUN = function(x) max(x, na.rm = TRUE),
            outer = ~0, inner = ~0,
            labels = list(), units = list())
convertFromS3groupedData(from)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ ~~Describe \code{formula} here~~ }
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{order.groups}{ ~~Describe \code{order.groups} here~~ }
  \item{FUN}{ ~~Describe \code{FUN} here~~ }
  \item{outer}{ ~~Describe \code{outer} here~~ }
  \item{inner}{ ~~Describe \code{inner} here~~ }
  \item{labels}{ ~~Describe \code{labels} here~~ }
  \item{units}{ ~~Describe \code{units} here~~ }
  \item{from}{ an object of class "groupedData" as used in package
    \code{nlme}
  }
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
\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{manip}


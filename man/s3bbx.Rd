\name{s3bbx}
\alias{s3bbx}
\encoding{latin1}
\non_function{}
\title{Covariates in the Rodriguez and Goldman simulation}
\usage{data(s3bbx)}
\description{
The \code{s3bbx} data frame has 2449 rows and 6 columns of the
covariates in the simulation by Rodriguez and Goldman of multilevel
dichotomous data.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{child}{a numeric vector identifying the child}
    \item{family}{a numeric vector identifying the family}
    \item{community}{a numeric vector identifying the community}
    \item{chldcov}{a numeric vector of the child-level covariate}
    \item{famcov}{a numeric vector of the family-level covariate}
    \item{commcov}{a numeric vector of the community-level covariate}
  }
}
\source{
  \url{http://data.princeton.edu/multilevel/simul.htm}
}
\references{
  Rodriguez, Germán and Goldman, Noreen (1995), "An assessment of
  estimation procedures for multilevel models with binary responses",
  \emph{Journal of the Royal Statistical Society, Series A}, \bold{158},
  73-89.
}
\examples{
data(s3bbx)
str(s3bbx)
}
\keyword{datasets}

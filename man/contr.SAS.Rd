\name{contr.SAS}
\alias{contr.SAS}
\title{Contrast Matrices for SAS-like Results}
\description{
  Return a matrix of contrasts suitable for comparison with SAS output.
}
\usage{
contr.SAS(n, contrasts = TRUE)
}
\arguments{
  \item{n}{a vector of levels for a factor, or the number of levels.}
  \item{contrasts}{a logical indicating whether contrasts should be
    computed.}
}
\details{
  This functions creates contrast matrices for use in
  fitting analysis of variance and regression models.  It is a trivial
  modification of \code{\link[stats]{contr.treatment}} that sets the
  base level to be the last level instead of the first level.  The
  coefficients produced when using these contrasts should be equivalent to
  those produced by SAS.
}
\value{
  A matrix with \code{n} rows and \code{k} columns, with \code{k=n-1} if
  \code{contrasts} is \code{TRUE} and \code{k=n} if \code{contrasts} is
  \code{FALSE}.
}
\references{
  Chambers, J. M. and Hastie, T. J. (1992)
  \emph{Statistical models.}
  Chapter 2 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth \& Brooks/Cole.
}
\seealso{
  \code{\link[stats]{contrasts}}, \code{\link[stats]{contr.treatment}}
}
\examples{
(cT <- contr.SAS(5))
all(crossprod(cT) == diag(4)) # TRUE: even orthonormal
}
\keyword{design}
\keyword{regression}
\keyword{array}

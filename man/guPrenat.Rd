\name{guPrenat}
\alias{guPrenat}
\docType{data}
\title{Prenatal care in Guatemala}
\description{
  Data on the prenatal care received by mothers in Guatemala.
}
\usage{data(guPrenat)}
\format{
  A data frame with 2449 observations on the following 15 variables.
  \describe{
    \item{kid}{a factor identifying the birth}
    \item{mom}{a factor identifying the mother or family}
    \item{cluster}{a factor identifying the community}
    \item{prenat}{a factor indicating if traditional or modern prenatal
      care was provided for the birth.}
    \item{childAge}{an ordered factor of the child's age at the time of
      the survey. }
    \item{motherAge}{a factor indicating if the mother was older or
      younger.  The cut-off age is 25 years.}
    \item{birthOrd}{an ordered factor for the birth's order within the family.}
    \item{indig}{a factor indicating if the mother is Ladino, or
      indigenous not speaking Spanish, or indigenous speaking Spanish.}
    \item{momEd}{a factor describing the mother's level of eduation.}
    \item{husEd}{a factor describing the husband's level of education.}
    \item{husEmpl}{a factor describing the husband's employment status.}
    \item{toilet}{a factor indicating if there is a modern toilet in the
      house.}
    \item{TV}{a factor indicating if there is a TV in the house and, if
      so, the frequency with which it is used.}
    \item{pcInd81}{the percentage of indigenous population in the
      community at the 1981 census.}
    \item{ssDist}{distance from the community to the nearest clinic.}
  }
}
%\details{}
\source{
  These data are available at
  \url{http://data.princeton.edu/multilevel/guPrenat.dat}.  Multiple
  indicator columns in the original data table have been collapsed to
  factors for this data frame.
}
\references{
  Rodriguez, Germán and Goldman, Noreen (1995), "Improved estimation
  procedures for multilevel models with binary response: a case-study",
  \emph{Journal of the Royal Statistical Society, Series A}, \bold{164},
  339-355.
}
\examples{
data(guPrenat)
summary(guPrenat)
}
\keyword{datasets}

\name{guImmun}
\encoding{latin1}
\alias{guImmun}
\docType{data}
\title{Immunization in Guatemala}
\description{
  Immunizations received by children in Guatemala.
}
\usage{data(guImmun)}
\format{
  A data frame with 2159 observations on the following 13 variables.
  \describe{
    \item{kid}{a factor identifying the child}
    \item{mom}{a factor identifying the family.}
    \item{comm}{a factor identifying the community.}
    \item{immun}{a factor indicating if the child received a
      complete set of immunizations. All children in this data frame
      received at least one immunization.}
    \item{kid2p}{a factor indicating if the child was 2 years or older
      at the time of the survey.}
    \item{mom25p}{a factor indicating if the mother was 25 years or older.}
    \item{ord}{an factor indicating the child's birth's order within the
      family.  Levels are \code{01} - first child, \code{23} - second or
      third child, \code{46} - fourth to sixth child, \code{7p} -
      seventh or later child.}
    \item{ethn}{a factor indicating the mother's ethnicity.  Levels are
      \code{L} - Ladino, \code{N} - indigenous not speaking Spanish, and
      \code{S} - indigenous speaking Spanish.}
    \item{momEd}{a factor describing the mother's level of eduation.
      Levels are \code{N} - not finished primary, \code{P} - finished
      primary, \code{S} - finished secondary}
    \item{husEd}{a factor describing the husband's level of education.
      Levels are the same as for \code{momEd} plus \code{U} - unknown.}
    \item{momWork}{a factor indicating if the mother had ever
      worked outside the home.}
    \item{rural}{a factor indicating if the family's location is
      considered rural or urban.}
    \item{pcInd81}{the percentage of indigenous population in the
      community at the 1981 census.}
  }
}
%\details{}
\source{
  These data are available at
  \url{http://data.princeton.edu/multilevel/guImmun.dat}.  Multiple
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
data(guImmun)
summary(guImmun)
}
\keyword{datasets}

\name{Early}
\alias{Early}
\docType{data}
\title{Early childhood intervention study}
\description{
  Cognitive scores of infants in a study of early childhood
  intervention.  The 103 infants from low income African American
  families were divided into a treatment group (58
  infants) and a control group (45 infants).  Starting at 0.5 years of age
  the infants in the treatment group were exposed to an enriched environment.
  Each infant's cognitive score on an age-specific, normalized scale was
  recorded at ages 1, 1.5, and 2 years.
}
\usage{data(Early)}
\format{
  This \code{groupedData} object contains the following columns
  \describe{
    \item{id}{An ordered factor of the id number for each infant.}
    \item{cog}{A numeric cognitive score.}
    \item{age}{The age of the infant at the measurement.}
    \item{trt}{A factor with two levels, \code{"N"} and \code{"Y"},
      indicating if the infant is in the early childhood intervention
      program.}
    }
}
%\source{}
\references{
  Singer, Judith D. and Willett, John B. (2003), \emph{Applied
    Longitudinal Data Analysis}, Oxford University Press. (Ch. 3)
}
\examples{
data(Early)
str(Early)
}
\keyword{datasets}

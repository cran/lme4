\name{egsingle}
\alias{egsingle}
\docType{data}
\title{US Sustaining Effects study}
\description{
  A subset of the mathematics scores from the U.S. Sustaining Effects
  Study.  The subset consists of information on 1721 students from 60 schools
}
\usage{data(egsingle)}
\format{
  A data frame with 7230 observations on the following 12 variables.
  \describe{
    \item{schoolid}{a factor of school identifiers}
    \item{childid}{a factor of student identifiers}
    \item{year}{a numeric vector indicating the year of the test}
    \item{grade}{a numeric vector indicating the student's grade}
    \item{math}{a numeric vector of test scores on the IRT scale score metric}
    \item{retained}{a factor with levels \code{0} \code{1} indicating if
    the student has been retained in a grade.}
    \item{female}{a factor with levels \code{Female} \code{Male}
      indicating the student's sex}
    \item{black}{a factor with levels \code{0} \code{1} indicating if
      the student is Black}
    \item{hispanic}{a factor with levels \code{0} \code{1} indicating if
    the student is Hispanic}
    \item{size}{a numeric vector indicating the number of students
      enrolled in the school}
    \item{lowinc}{a numeric vector giving the percentage of low-income
      students in the school}
    \item{mobility}{a numeric vector}
  }
}
%\details{}
\source{
  These data are distributed with the HLM software package (Bryk,
  Raudenbush and Congdon, 1996).  Conversion to the R format is
  described in Doran and Lockwood (2004).
}
\references{
  Doran, Harold C. and Lockwood, J.R. (2004), \emph{Fitting value-added
    models in R}, (submitted).
}
\examples{
data(egsingle)
str(egsingle)
fm1 <- lme(math~year*size+female, egsingle,
           random = ~1|schoolid/childid,
           control = list(EMv = TRUE, msV = TRUE, opt = 'optim'))
summary(fm1)
}
\keyword{datasets}

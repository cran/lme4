\name{Hsb82}
\alias{Hsb82}
\docType{data}
\title{High School and Beyond - 1982}
\description{
  Data from the 1982 study ``High School and Beyond''.
}
\usage{data(Hsb82)}
\format{
  A data frame with 7185 observations on students including the
  following 8 variables.
  \describe{
    \item{school}{an ordered factor designating the school that the
      student attends.}
    \item{minrty}{a factor with levels}
    \item{sx}{a factor with levels \code{Male} and \code{Female}}
    \item{ses}{a numeric vector of socio-economic scores}
    \item{mAch}{a numeric vector of Mathematics achievement scores}
    \item{meanses}{a numeric vector of mean \code{ses} for the school}
    \item{sector}{a factor with levels \code{Public} and \code{Catholic}}
    \item{cses}{a numeric vector of centered \code{ses} values where the
    centering is with respect to the \code{meanses} for the school.}
  }
}
\details{
  Each row in this data frame contains the data for one student.
}
%\source{}
\references{
  Raudenbush, Stephen and Bryk, Anthony (2002), \emph{Hierarchical
    Linear Models: Applications and Data Analysis Methods}, Sage
  (chapter 4).
}
\examples{
data(Hsb82)
summary(Hsb82)
}
\keyword{datasets}

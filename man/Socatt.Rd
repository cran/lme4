\name{Socatt}
\alias{Socatt}
\docType{data}
\title{Social Attitudes Survey}
\description{
  These data come from the British Social Attitudes (BSA) Survey started
  in 1983. The eligible persons were all adults aged 18 or over living
  in private households in Britain.  The data consist of completed
  results of 264 respondents out of 410.  
}
\usage{data(Socatt)}
\format{
  A data frame with 1056 observations on the following 9 variables.
  \describe{
    \item{district}{District ID - a factor}
    \item{respond}{Respondent code (within district) - a factor}
    \item{year}{A factor with levels \code{1983}, \code{1984},
      \code{1985}, and \code{1986}}
    \item{numpos}{An ordered factor giving the number of positive answers to
      seven questions.}
    \item{party}{Political party chosen - a factor.  Levels are
      \code{conservative}, \code{labour}, \code{Lib/SDP/Alliance},
      \code{others}, and \code{none}.}
    \item{class}{Self assessed social class - a factor.  Levels are
      \code{middle}, \code{upper working}, and \code{lower working}.}
    \item{gender}{Respondent's sex. (1=male, 2=female)}
    \item{age}{Age in years}
    \item{religion}{Religion - a factor.  Levels are \code{Roman
	Catholic}, \code{Protestant/Church of England}, \code{others},
       and \code{none}.}
  }
}
\details{
  These data are provided at \url{http://multilevel.ioe.ac.uk/softrev/}
  as an example of multilevel data with a multinomial response.  Results
  from several multilevel software packages are provided at this site.
}
\source{
  \url{http://multilevel.ioe.ac.uk/softrev/socatt.htm}
}
\references{
  McGrath, K. and Waterton, J. (1986). \emph{British Social Attitudes
    1983-1986 panel survey.} London, Social and Community Planning
  Research.
}
\examples{
data(Socatt)
str(Socatt)
summary(Socatt)
}
\keyword{datasets}

\name{groupedData-class}
\docType{class}
\alias{groupedData-class}
\alias{$,groupedData-method}
\alias{$<-,groupedData-method}
\alias{[,groupedData,missing,ANY,logical-method}
\alias{[,groupedData,missing,ANY,missing-method}
\alias{[,groupedData,ANY,missing,logical-method}
\alias{[,groupedData,ANY,missing,missing-method}
\alias{[,groupedData,ANY,ANY,missing-method}
\alias{[,groupedData,ANY,ANY,logical-method}
\alias{[<-,groupedData,ANY,ANY-method}
\alias{[<-,groupedData,ANY,missing-method}
\alias{[<-,groupedData,missing-method}
\alias{[[,groupedData,ANY,missing-method}
\alias{[[,groupedData,ANY,ANY-method}
\alias{coerce,groupedData,data.frame-method}
\alias{formula,groupedData-method}
\alias{lmList,formula,groupedData-method}
\alias{lmList,missing,groupedData-method}
\alias{plot,groupedData,missing-method}
\alias{show,groupedData-method}
\alias{summary,groupedData-method}

\title{Class "groupedData"}
\description{A \code{groupedData} object is at least a data.frame and a
  display  formula that provides a grouping structure.  At present we
  are using an old-style class for this.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("groupedData",
    ...)} or by calls to the \code{groupedData} constructor function.
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{The actual data; an object of class \code{"data.frame"} }
    \item{\code{formula}:}{The display formula; an object of class
      \code{"formula"}} 
    \item{\code{outer}:}{A formula all of whose variables are outer to the
      grouping factor}
    \item{\code{inner}:}{A formula of factors all of whose variables are
      inner to the grouping factor.}
    \item{\code{labels}:}{An optional list of default labels for the
      \code{x} and \code{y} axes of the default plot.}
    \item{\code{units}:}{An optional list of default units
      specifications for the \code{x} and \code{y} axes of the default
      plot.} 
  }
}
\section{Methods}{
  \describe{
    \item{\$}{\code{signature(x = "groupedData")}: ... }
    \item{\$<-}{\code{signature(x = "groupedData")}: ... }
    \item{[}{\code{signature(x = "groupedData", i = "missing", j = "ANY", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "groupedData", i = "missing", j = "ANY", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "groupedData", i = "ANY", j = "missing", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "groupedData", i = "ANY", j = "missing", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "groupedData", i = "ANY", j = "ANY", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "groupedData", i = "ANY", j = "ANY", drop = "logical")}: ... }
    \item{[<-}{\code{signature(x = "groupedData", i = "ANY", j = "ANY")}: ... }
    \item{[<-}{\code{signature(x = "groupedData", i = "ANY", j = "missing")}: ... }
    \item{[<-}{\code{signature(x = "groupedData", i = "missing")}: ... }
    \item{[[}{\code{signature(x = "groupedData", i = "ANY", j = "missing")}: ... }
    \item{[[}{\code{signature(x = "groupedData", i = "ANY", j = "ANY")}: ... }
    \item{coerce}{\code{signature(from = "groupedData", to = "data.frame")}: ... }
    \item{formula}{\code{signature(x = "groupedData")}: ... }
    \item{lmList}{\code{signature(formula = "formula", data = "groupedData")}: ... }
    \item{lmList}{\code{signature(formula = "missing", data = "groupedData")}: ... }
    \item{plot}{\code{signature(x = "groupedData", y = "missing")}: ... }
    \item{show}{\code{signature(object = "groupedData")}: ... }
    \item{summary}{\code{signature(object = "groupedData")}: ... }
  }
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link{groupedData}}}
%\examples{}
\keyword{classes}

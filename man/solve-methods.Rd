\name{solve-methods}
\docType{methods}
\alias{solve-methods}
\alias{solve,pdBlocked,missing-method}
\alias{solve,pdDiag,missing-method}
\alias{solve,pdIdent,missing-method}
\alias{solve,pdLogChol,missing-method}
\alias{solve,pdMat,missing-method}
\alias{solve,pdNatural,missing-method}
%\alias{solve,pdScalar,missing-method}
\title{Methods for Function solve}
\description{The \code{solve} methods for \code{pdMat} classes all use
  a single argument call.  They return a matrix of the same class as
  \code{a}.
  }
\section{Methods}{
  \describe{
    \item{a = pdBlocked, b = missing}{invert a \code{pdBlocked} object}
    \item{a = pdDiag, b = missing}{invert a \code{pdDiag} object}
    \item{a = pdIdent, b = missing}{invert a \code{pdIdent} object}
    \item{a = pdLogChol, b = missing}{invert a \code{pdLogChol} object}
    \item{a = pdMat, b = missing}{invert a \code{pdMat} object}
    \item{a = pdNatural, b = missing}{invert a \code{pdNatural} object}
  }
}
\keyword{manip}
\keyword{algebra}
\keyword{array}
\keyword{methods}

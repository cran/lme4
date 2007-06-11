\name{nlmer-class}
\docType{class}
\alias{nlmer-class}
\alias{show,nlmer-method}
%\alias{VarCorr,nlmer-method}

\title{Representation of a Nonlinear Mixed Model} 
\description{The \code{lmer} class is a representaiton of nonlinear
  mixed model using sparse matrices.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("nlmer", ...)}
  or, more commonly via the \code{\link{nlmer}} function.%
}
\section{Slots}{
  \describe{
    \item{\code{env}:}{An environment (of class \code{"environment"})
      created for the evaluation of the nonlinear model function.}
    \item{\code{model}:}{The nonlinear model function as an object of
      class \code{"call"}.}
    \item{\code{frame}:}{The model frame for the model, as an object of
      class \code{"data.frame"}.}
    \item{\code{pnames}:}{Names of the parameters in the nonlinear model
      (class \code{"character"}).}
    \item{\code{call}:}{The matched call to the function that
      created the object. (class \code{"call"}).} 
    \item{\code{flist}:}{The list of grouping factors for the random
      effects.}
    \item{\code{Xt}:}{Sparse model matrix (class \code{"dgCMatrix"}) for
      the fixed effects.}
    \item{\code{Zt}:}{Sparse model matrix (class \code{"dgCMatrix"}) for
      the random effects.}
    \item{\code{y}:}{The response vector (class \code{"numeric"}).}
    \item{\code{weights}:}{Numeric weights vector. This may be of length
      zero (0) indicating unit weights.}
    \item{\code{cnames}:}{a list of character vectors of column names
      of the random-effects variance-covariance matrix associated with
      each grouping factor and the fixed-effects model matrix.}
    \item{\code{Gp}:}{integer vector of group pointers within the random
      effects and the rows of the transposed model matrix in the
      \code{Zt} slot.  The elements of \code{Gp} are the 0-based index
      of the first element corresponding to each grouping factor.  Thus
      the first element is always 0.  The last element is the total
      length of the random effects vector (also the total number of rows
      in the matrix in the \code{Zt} slot).}
    \item{\code{dims}:}{A named integer vector of dimensions.}
    \item{\code{ST}:}{A list of S and T factors in the TSST' Cholesky
      factorization of the relative variance matrices of the random
      effects associated with each grouping factor.  The unit lower
      triangular matrix T and the diagonal matrix S are stored as a
      single matrix whose diagonal elements determine S while the
      subdiagonal elements determine the non-trivial elements of T.}
    \item{\code{Vt}:}{A sparse matrix (of class  \code{"dgCMatrix"})
      containing the transpose of V=ZTS.}
    \item{\code{L}:}{The Cholesky decomposition (class
      \code{"CHMfactor"}) of V'V+I.}
    \item{\code{mu}:}{A numeric vector of predicted values from the
      model with an attibute \code{"gradient"}, the gradient matrix.}
    \item{\code{Mt}:}{The derivative of \code{mu} with respect to
      \code{uvec} stored as a sparse matrix (class \code{"dgCMatrix"}).}
    \item{\code{deviance}:}{Named numeric vector of containing the
      deviance corresponding to the maximum likelihood (ML) and REML
      criteria and various components.}
     \item{\code{fixef}:}{Numeric vector of fixed effects.}
    \item{\code{ranef}:}{Numeric vector of random effects on the
      original scale.}
    \item{\code{uvec}:}{Numeric vector of orthogonal, constant variance
      random effects.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "nlmer")}}
    \item{VarCorr}{\code{signature(x = "nlmer")}: Extract the variances,
      standard deviations and correlations of the random effects.}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link{nlmer}}, \code{\linkS4class{lmer}}}
\examples{
showClass("nlmer")
}
\keyword{classes}

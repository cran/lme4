### pdLogChol - a general positive definite structure parameterized
###   by the non-zero elements of the Cholesky factor.  The logarithms of
###   the diagonal elements are the first Ncol elements of the parameter
###   vector

setGeneric('pdLogChol',
           function(value, form=formula(NULL), nam = character(), data=list(),
                    ...)
           standardGeneric('pdLogChol'))

setMethod('pdLogChol',
          signature(value = 'formula', form = 'missing',
                    nam = 'missing', data = 'missing'),
          function(value, form, nam, data, ...) {
              new('pdLogChol', form = value)
          })

## Methods for the pdLogChol class

setReplaceMethod("coef",
                 signature(object = "pdLogChol", value = "numeric"),
                 function(object, value) {
                     .Call("pdLogChol_coefGets", object, value, PACKAGE = "lme4")
#                     npar <- length(value)
#                     if (npar < 1)
#                         stop('coef for a pdLogChol object must have length > 0')
#                     if (npar != length(object@param)) {
#                         Ncol <- round((sqrt(8*length(value) + 1) - 1)/2)
#                         np <- (Ncol * (Ncol + 1))/2
#                         if (np != npar)
#                             stop(paste("coef for a pdLogChol object cannot have",
#                                        "length", npar))
#                         lenPar <- length(object@param)
#                         if (lenPar <= 0 && length(object@Names) > 0) {
#                             lenPar <- length(object@Names)
#                             lenPar <- (lenPar * (lenPar+1))/2
#                         }
#                         if (lenPar && lenPar != npar)
#                             stop("coef for a pdLogChol object has inconsistent length")
#                         object@Ncol <- as.integer(Ncol)
#                         object@factor <- matrix(0., Ncol, Ncol)
#                     }
#                     Ncol <- object@Ncol
#                     fact <- object@factor
#                     diag(fact) <- exp(value[1:Ncol])
#                     fact[row(fact) < col(fact)] <- value[-(1:Ncol)]
#                     object@param <- value
#                     object@factor <- fact
#                     object@logDet <- sum(value[1:Ncol])
#                     object
                 })

setAs("pdLogChol", "pdmatrix",
      function(from) new("pdmatrix", crossprod(from@factor)),
      function(from, value) {
          nc <- ncol(value)
          if (!identical(nc, dim(value)[2]))
              stop("value must be a square matrix")
          if (length(from@param) < 1) {
              from@Ncol <- nc
          }
          if (from@Ncol != nc)
              stop("can not change length of an initialized pdMat object")
          Names <- dimnames(value)[[2]]
          if (!is.null(Names))
              from@Names <- Names
          fact <- .Call("nlme_Chol", as(value, "matrix"), PACKAGE="lme4")
          from@factor <- fact
          from@logDet = sum(log(diag(fact)))
          from@param <- c(log(diag(fact)), fact[col(fact) > row(fact)])
          from
      })

setMethod("solve", signature(a="pdLogChol", b="missing"),
          function(a, b) {
              if (!isInitialized(a))
                  stop(paste("Uninitialized", class(a), "object"))
              as(a, "pdmatrix") <- crossprod(t(solve(a@factor)))
              a
          })

#setMethod("summary", signature(object="pdLogChol"),
#          function(object, structName, noCorrelation, ...) {
#              if (missing(structName)) structName =
#                   "General positive-definite, Log-Cholesky parametrization"
#              if (missing(noCorrelation)) noCorrelation = FALSE
#              callNextMethod()
#          })

setMethod("LMEgradient",
          signature(x="pdLogChol", A="matrix", nlev="numeric"),
          function(x, A, nlev)
          .Call("pdLogChol_LMEgradient", x, A, nlev, PACKAGE="lme4")
          )

setMethod("LMEhessian",
          signature(x="pdLogChol", A="matrix", H="array",
                    nlev="numeric"),
          function(x, A, H, nlev)
      {
          .Call("pdLogChol_LMEhessian", x, A, H, nlev, PACKAGE="lme4")
      })

setReplaceMethod("EMupdate",
                 signature(x="pdLogChol", nlev="numeric", value="matrix"),
                 function(x, nlev, value) {
                     .Call("pdLogChol_EMupdate", x, nlev, value, PACKAGE="lme4")
#                     if (!isInitialized(x))
#                         stop(paste("Uninitialized", class(x), "object"))
#                     if (any(dim(value) != dim(x)))
#                         stop(paste("value must be a matrix of dimension",
#                                    dim(x)))
#                     if (length(nlev) != 1 || nlev <= 0)
#                         stop("nlev must be > 0")
#                     as(x, "pdmatrix") <- nlev*crossprod(t(solve(value)))
#                     x
                 })

setMethod("pdgradient", "pdLogChol",
          function(x) {
              .Call("pdLogChol_pdgradient", x, PACKAGE="lme4")
#              fact <- as(x, "pdfactor")
#              pars <- x@param
#              Ncol <- ncol(fact)
#              dn <- dimnames(fact)
#              if (!is.null(dn)) dn <- c(list(dimnames(fact)), NULL)
#              val <- array(0., dim = c(dim(fact), length(pars)),
#                           dimnames = dn)
#              nc <- 1:Ncol
#              val[cbind(nc, nc, nc)] <- exp(pars[nc])
#              offdiag <- row(fact) < col(fact)
#              val[cbind(row(fact)[offdiag], col(fact)[offdiag],
#                  (Ncol+1):length(pars))] <- 1
#              for (i in seq(along = pars)) {
#                  pr <- crossprod(fact, val[,,i])
#                  val[,,i] <- pr + t(pr)
#              }
#              val
          })


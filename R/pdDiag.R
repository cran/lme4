### pdDiag - diagonal structure parameterized by the logarithm of
###   the square root of the diagonal terms.

setGeneric('pdDiag',
           function(value, form, nam, data, ...)
           standardGeneric('pdDiag'))

setMethod("pdDiag",
          signature(value = 'formula', form = 'missing',
                    nam = 'missing', data = 'missing'),
          function(value, form, nam, data, ...) {
              new('pdDiag', form = value)
          })

## Methods for the pdDiag class

setReplaceMethod("coef",
                 signature(object = "pdDiag", value = "numeric"),
                 function(object, value) {
                     .Call("pdDiag_coefGets", object, value, PACKAGE = "lme4")
#                     lenVal <- length(value)
#                     if (lenVal <= 0)
#                         stop('coef for a pdDiag object must have length > 0')
#                     lenPar <- length(object@param)
#                     if (lenPar == 0)
#                         lenPar <- length(object@Names)
#                     if (lenPar != 0 && lenPar != lenVal)
#                         stop("coef for a pdDiag object must be same as its number of rows")
#                     object@param <- value
#                     object@Ncol <- lenVal
#                     object@factor <- diag(exp(value), ncol = lenVal)
#                     object@logDet <- sum(value)
#                     object
                 })

setAs("pdDiag", "pdmatrix",
      function(from) {
          if (!isInitialized(from))
              stop(paste("Uninitialized", class(from), "object"))
          value <- diag(exp(2 * from@param), ncol = from@Ncol)
          nam <- from@Names
          if (length(nam) == length(from@param)) {
              dimnames(value) <- list(nam, nam)
          }
          new("pdmatrix", value)
      },
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
          coef(from) <- log(diag(value))/2
          from
      })

setMethod("solve", signature(a="pdDiag", b="missing"),
          function(a, b) {
              if (!isInitialized(a))
                  stop(paste("Uninitialized", class(a), "object"))
              coef(a) <- -a@param
              a
          })

setMethod("LMEgradient",
          signature(x="pdDiag", A="matrix", nlev="numeric"),
          function(x, A, nlev) {
              .Call("pdDiag_LMEgradient", x, A, nlev, PACKAGE="lme4")
          })

setMethod("LMEhessian",
          signature(x="pdDiag", A="matrix", H="array",
                    nlev="numeric"),
          function(x, A, H, nlev)
      {
          theta <- x@param
          q <- length(theta)

          ## the part involving D_iD_j
          ans <-
              if (q > 1)
                  diag(-exp(2*theta)*colSums(A*A))
              else as.matrix(-exp(2*theta)*A*A)

          ## add the part not involving D_iD_j
          tmp <- matrix(0, nrow=q, ncol=q)
          for (j in 1:q)
              tmp[, j] <- H[seq(from = 1+(j-1)*q*q*(q+1), length=q, by = q+1)]
          ans <- ans + tmp*outer(theta, theta,
                                 function(v, w) exp(2*(v+w)))
          nm <- names(x)
          ans <- ans+t(ans)
          if (!is.null(nm))
              dimnames(ans) <- list(nm, nm)
          ans
      })

setReplaceMethod("EMupdate",
                 signature(x="pdDiag", nlev="numeric", value="matrix"),
                 function(x, nlev, value) {
                     .Call("pdDiag_EMupdate", x, nlev, value, PACKAGE="lme4")
                 })

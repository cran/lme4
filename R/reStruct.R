###      Methods for the class of random-effects structures.

setMethod("reStruct", signature(fixed = "formula",
                                random = "list",
                                data = "data.frame",
                                weights = "missing",
                                REML = "logical"),
          function(fixed, random, data, weights, REML, nextraCols)
      {
          reStruct(fixed, random, data, weights=numeric(0), REML, nextraCols)
      })

setMethod("reStruct", signature(fixed = "formula",
                                random = "list",
                                data = "data.frame",
                                weights = "numeric",
                                REML = "logical"),
          function(fixed, random, data, weights, REML, nextraCols)
      {
          ## given a matrix of no. of rows needed in the last row of a
          ## level for each column, create the storedRows or
          ## decomposedRows list
          createRows <-
              function(rowMatrixEnd)
              {
                  rowMatrixBeg <- as.vector(rowMatrixEnd)
                  nonzero <- rowMatrixBeg[rowMatrixBeg != 0]
                  rowMatrixBeg[rowMatrixBeg != 0] <-
                      c(1, nonzero[-length(nonzero)]+1)
                  dim(rowMatrixBeg) <- dim(rowMatrixEnd)
                  value <- lapply(seq(length=nrow(rowMatrixBeg)),
                                  function(i)
                              {
                                  ind <- rowMatrixBeg[i, ] != 0
                                  beg <- rowMatrixBeg[i, ind]
                                  end <- rowMatrixEnd[i, ind]
                                  lapply(seq(along=beg),
                                         function(i) seq(from=beg[i],
                                                         to=end[i]))
                              })
                  names(value) <- dimnames(rowMatrixEnd)[[1]]
                  value
              }
          ## order from inner to outer groups
          random <- rev(random)
          self <- new("reStruct",
                      fixed=fixed,
                      REML=REML)

          ## Get the grouping factors
          self@groups <- as.data.frame(lapply(names(random),
                                              function(expr)
                                              eval(parse(text=expr),
                                                   data,
                                                   environment(fixed))))
          for (i in seq(to=1, length=length(self@groups)-1, by = -1)) {
              self@groups[, i] <-
                  as.factor(paste(as.character(self@groups[[i+1]]),
                                  as.character(self@groups[[i]]),
                                  sep = "/"))
          }
          names(self@groups) <- names(random)
          ## save the old ordering
          row.names(self@groups) <- row.names(data)
          ## generate the new ordering for the model.frame
          self@origOrder <- do.call("order", self@groups)
          self@groups[, 1] <- as.factor(self@groups[,1])
          self@reverseOrder <- order(self@origOrder)
          ## Reorder the model.frame
          data[,] <- data[self@origOrder, , drop = FALSE]
          row.names(data) <- row.names(self@groups)[self@origOrder]
          self@groups[,] <- self@groups[self@origOrder, , drop = FALSE]
          row.names(self@groups) <- row.names(data)
          ## Create the model.matrix for fixed and response
          self@original <- model.matrix(fixed, data)
          self@assign.X <- attr(self@original, "assign")
          self@original <- cbind(self@original, model.response(data))
          nCol <- ncol(self@original)
          p <- nCol-1
          ## Set the column name for response
          colnames(self@original)[nCol] <-
              colnames(data)[attr(attr(data, "terms"),
                                        "response")]

          qVector <- integer(length(random))
          for (i in seq(length=length(random),
                        by = -1, to = 1)) {
              self@original <- cbind(model.matrix(formula(random[[i]]),
                                                  data),
                                     self@original)
              qVector[i] <- ncol(self@original)-nCol
              nCol <- ncol(self@original)
          }

          nextraCols <- as.integer(nextraCols)
          if (nextraCols > 0) {
              self@original <- cbind(self@original,
                                     matrix(as.numeric(NA),
                                            nrow=nrow(self@original),
                                            ncol=nextraCols))
              nCol <- ncol(self@original)
          }

          ncols <- qVector
          qVector <- c(qVector, p, 1)
          indx <- seq(along=random)
          names(indx) <- names(random)
          reStructColumns <- lapply(seq(along=qVector),
                                    function(i, end=cumsum(qVector))
                                    as.integer(seq(length=qVector[i],
                                                   to=end[i])))
          if (length(weights) > 0) {
              weighted(self) <- weights[self@origOrder]
              rm(weights)
          }
          Q <- ncol(self@groups)
          random <- lapply(indx, function(i)
                           lmeLevel(random[[i]],
                                    self@groups[[i]],
                                    reStructColumns[[i]],
                                    if (self@useWeighted)
                                    self@weighted
                                    else self@original))
          class(random) <- "lmeLevelList"
          self@random <- random
          rm(random)
          ## FIXME - there should be a less grubby way of doing this
          pars <- lapply(self@random, coef)
          maxInd <- 0
          for (i in seq(along = pars)) {
              self@random[[i]]@parsInd <- as.integer(maxInd + seq(along = pars[[i]]))
              maxInd <- maxInd + length(pars[[i]])
          }

          self@random[["*fixed*"]] <- lmeLevel(columns=reStructColumns[[Q+1]],
                                               modelMatrix=self@original)
          self@random[["*response*"]] <- lmeLevel(columns=reStructColumns[[Q+2]],
                                                  modelMatrix=self@original)
          N <- nrow(self@groups)

          self <- .Call("nlme_reStructDims", self, PACKAGE="lme4")

          self@decomposed <- self@stored <-
              matrix(0, nrow = self@random[[Q+2]]@storedRows[[1]], ncol=nCol)
          self@bbetas <- numeric(nrow(self@stored))

          self
      })

setMethod("coef", "reStruct",
          function(object, ...)
          unlist(lapply(object@random, coef)))

setReplaceMethod("coef", signature(object="reStruct",
                                   value="numeric"),
                 function(object, value)
             {
                 names(value) = NULL
                 for (i in seq(length=length(object@random)-2)) {
                     coef(object@random[[i]]) =
                         value[object@random[[i]]@parsInd]
                 }
                 object@logLik = as.numeric(NA)
                 object
             })

setMethod("getResponse", signature(object="reStruct"),
          function(object, form)
      {
          object@original[, object@random[["*response*"]]@columns]
      })

setReplaceMethod("response", signature(x="reStruct", value="numeric"),
                 function(x, value)
             {
                 if (length(value) != nrow(x@original))
                     stop("Dimension mismatch in model.matrix")
                 x@original[, x@random[["*response*"]]@columns] = value
                 x@dirtyDecomposed = TRUE
                 x@useWeighted = FALSE
                 x
             })

setReplaceMethod("weighted", signature(x="reStruct",
                                       value="numeric"),
                 function(x, value)
             {
                 if (all(value == 1.0)) {
                     if (x@useWeighted) {
                         x@useWeighted = FALSE
                         x@dirtyDecomposed = TRUE
                     }
                 } else {
                     ## make sure the row number of original
                     ## model.matrix match length of weights
                     if (nrow(x@original) != length(value) &&
                         length(value) != 1)
                         stop("Dimension mismatch in setting weighted")
                     if (identical(dim(x@weighted), x@original))
                         x@weighted[] = x@original * value
                     else x@weighted = x@original * value
                     x@useWeighted = TRUE
                     x@dirtyDecomposed = TRUE
                 }
                 x
             })
setMethod("model.matrix", "reStruct",
          function(object, ...) object@original)

setReplaceMethod("model.matrix", signature(x="reStruct",
                                           value="matrix"),
                 function(x, value)
             {
                 ## FIXME: check that value is a model.matrix

                 ## make sure the dimensions of old and new
                 ## model.matrix match
                 if (!identical(dim(value), dim(x@original)))
                     stop("Dimension mismatch in model.matrix")
                 x@original = value
                 x@dirtyDecomposed = TRUE
                 x@useWeighted = FALSE
                 x
             })

setMethod("getGroups", signature(object="reStruct",
                                 form="missing",
                                 data="missing",
                                 sep="missing"),
          function(object, form, level, data, sep)
      {
          Q <- length(object@random)-2
          if (missing(level))
              level <- Q
          else if (length(level) > Q || max(level) > Q)
              stop ("Invalid value for level")
          val <- lapply(object@random[seq(length=Q)[level]],
                        function(ranlev) ranlev@groups)
          if (length(val) == 1) { # single group
              valname <- names(val)
              val <- val[[1]]
              attr(val, "label") <- valname
              val
          } else {
              as.data.frame(val)
          }
      })

setMethod("logLik", signature(object="reStruct"),
          function(object)
      {
          value = .Call("nlme_logLikelihood", object,
                         NULL,          # do not pass new parameter value
                         PACKAGE="lme4")
          p = length(object@random[["*fixed*"]]@columns)
          # df calculated from sigma + fixed effects + random effects pars
          attr(value, "df") = 1 + p +
              sum(unlist(lapply(object@random, function(x)length(coef(x)))))
          attr(value, "nall") = nrow(object@original)
          attr(value, "nobs") = nrow(object@original) - p
          class(value) = "logLik"
          value
      })

setReplaceMethod("EMsteps", signature(x="reStruct", value="list"),
                 function(x, value)
             {
                 x <- .Call("nlme_reStructEMsteps", x, value$niterEM,
                            value$EMverbose,
                            PACKAGE="lme4")
#                  verbose = value$EMverbose
#                  randIndx = seq(length=length(x@random)-2)
#                  for (i in seq(length=value$niterEM)) {
#                      x = .Call("nlme_commonDecompose", x, NULL,
#                                 PACKAGE="lme4")
#                      for (j in randIndx)
#                          EMupdate(x@random[[j]]) =
#                              x@random[[j]]@updateFactor
#                      if (verbose) {
#                          cat("\n**EM Iteration", i,
#                              .Call("nlme_logLikelihood", x,
#                                    NULL, # do not pass new parameter value
#                                    PACKAGE="lme4"), coef(x), "\n")
#                      }
#                      x@logLik = as.numeric(NA)
#                  }
                 .Call("nlme_commonDecompose", x, NULL,
                       PACKAGE="lme4")
             })

setReplaceMethod("LMEoptimize", signature(x="reStruct",
                                          value="list"),
                 function(x, value)
             {
                 if (value$msMaxIter < 1)
                     return(x)
                 xval = -.Call("nlme_logLikelihood",
                                x,
                                NULL,  # use already set parameters
                                PACKAGE="lme4")
                 xval =
                     if (xval > 0)
                         xval+1
                     else abs(min(xval/2, xval+1))
                 if (value$optimizer == "nlm") {
                     optimRes =
                         if (value$analyticGradient) {
                             optim(fn = function(params)
                                   .Call("nlme_logLikelihood",
                                         x,
                                         params,  # new parameter value
                                         PACKAGE="lme4"),
                                   gr = function(params)
                                   LMEgradient(.Call("nlme_commonDecompose",
                                                     x, params,
                                                     PACKAGE="lme4")),
                                   par = c(coef(x)), #hessian = TRUE,
                                   method = "BFGS",
                                   control = list(trace = value$msVerbose,
                                   reltol = value$msTol,
                                   fnscale = -1,
#                                   fnscale = -xval,
#                                   parscale = 1/value$msScale(coef(x)),
                                   maxit = value$msMaxIter))
                         } else {
                             optim(fn = function(params)
                                   .Call("nlme_logLikelihood",
                                         x,
                                         params,  # new parameter value
                                         PACKAGE="lme4"),
                                   par = c(coef(x)), #hessian = TRUE,
                                   method = "BFGS",
                                   control = list(trace = value$msVerbose,
                                   reltol = value$msTol,
                                   fnscale = -1,
#                                   fnscale = -xval,
#                                   parscale = 1/value$msScale(coef(x)),
                                   maxit = value$msMaxIter))
                         }
                     if (optimRes$convergence != 0) {
                         warning("optim failed to converge")
                     }
                     .Call("nlme_commonDecompose", x, optimRes$par,
                           PACKAGE="lme4")
#                  } else if (value$optimizer == "ms") {
#                      pars <- coef(x)
#                      .Call("nlme_msOptimize", value$msMaxIter,
#                            value$msTol, rep(1.0, length(pars)),
#                            value$msVerbose, x, pars,
#                            value$analyticGradient,
#                            PACKAGE = "lme4")
                 } else {
#                     typsize <- 1/value$msScale(coef(x))
                     typsize <- rep(1.0, length(coef(x)))
                     if (is.null(value$nlmStepMax))
                         value$nlmStepMax <-
                             max(100 * sqrt(sum((coef(x)/typsize)^2)), 100)
                     nlmRes =
                         nlm(f = if (value$analyticGradient) {
                             function(params)
                             {
                                 x = .Call("nlme_commonDecompose",
                                            x, params,
                                            PACKAGE="lme4")
                                 grad = -LMEgradient(x)
                                 ans = -x@logLik
                                 attr(ans, "gradient") = grad
                                 ans
                             }
                         } else {
                             function(params)
                             {
                                 -.Call("nlme_logLikelihood",
                                        x,
                                        params,
                                        PACKAGE="lme4")
                             }
                         },
                             p = c(coef(x)), #hessian = TRUE,
                             print.level = if (value$msVerbose) 2 else 0,
                             steptol = value$msTol,
                             gradtol = value$msTol,
                             stepmax = value$nlmStepMax,
                             typsize=typsize,
#                             fscale=xval,
                             iterlim = value$msMaxIter)
                     .Call("nlme_commonDecompose", x, nlmRes$estimate,
                           PACKAGE="lme4")
                 }
             })


setMethod("LMEgradient", signature(x="reStruct", A="missing", nlev="missing"),
          function(x, A, nlev)
      {
          unlist(lapply(x@random[seq(length=length(x@random)-2)],
                        LMEgradient))
      })

setMethod("fitted", signature(object="reStruct"),
          function(object, ...)
      {
          .Call("nlme_reStruct_fitted", object, PACKAGE="lme4")
      })

setMethod("fixef", signature(object="reStruct"),
          function(object, ...) {
              fixd = object@random[['*fixed*']]
              val = object@bbetas[fixd@storedRows[[1]]]
              nn = dimnames(object@original)[[2]][fixd@columns]
              if (length(nn) == length(val)) names(val) = nn
              val
          })

setMethod("ranef", signature(object="reStruct"),
          function(object, ...) {
              lapply(object@random[-(length(object@random) - c(1,0))],
                     function(x) matrix(object@bbetas[unlist(x@storedRows)],
                                        ncol = length(x@columns),
                                        byrow = TRUE))
          })

setMethod("summary", "reStruct",
          function(object, ...) {
              fixd = object@random[['*fixed*']]
              fstrRows = fixd@storedRows[[1]]
              fcols = fixd@columns
              rsp = object@random[['*response*']]
              sigma = abs(object@bbetas[rsp@storedRows[[1]]])
              nobs = dim(object@original)[1]
              nfixed = length(fixd@columns)
              denomDF = nobs - ifelse(object@REML, nfixed, 0)
              sigma = sigma/sqrt(denomDF)
              fcoef = object@bbetas[fstrRows]
              DF = getFixDF(object)$X
              rinv =
                  .Call("nlme_commonDecompose", object, NULL,
                        PACKAGE="lme4")@stored[fstrRows, fcols, drop=FALSE]
              se = sqrt(rowSums(rinv * rinv))
              corF = new("corrmatrix", crossprod(t(rinv/se)), stdDev = se)
              coefs = cbind(fcoef, se, DF)
              nn = dimnames(object@original)[[2]][fcols]
              dimnames(coefs) =
                  list(nn, c("Estimate", "Std. Error", "DF"))
              rnd = rev(object@random)[-(1:2)]
              ngrps = sapply(rnd, function(x) x@nlev)
              new("summary.reStruct",
                  fixed = object@fixed,
                  coefficients = as.matrix(coefs),
                  scale = sigma,
                  denomDF = as.integer(denomDF),
                  REML = object@REML,
                  ngrps = ngrps,
                  nobs = nobs,
                  corFixed = corF,
                  reSumry = lapply(rnd,
                      function(x) summary(solve(x@precision))))
          })

setMethod("getFixDF", signature(object="reStruct"),
          function(object)
      {
          ## calculates degrees of freedom for fixed effects Wald tests
#          Q <- length(object@random)-2
#          columns = object@random[["*fixed*"]]@columns
#          X = object@original[, columns, drop = FALSE]
#          ngrps = unlist(lapply(object@random, function(lmeLevel)
#                                 lmeLevel@nlev))
#          names(ngrps) = names(object@random)
          val = .Call("nlme_getFixDF", object)
          names(val$X) =
              colnames(object@original)[object@random[["*fixed*"]]@columns]
                                        # Convert R's assign to S-PLUS style
          assign = object@assign.X
          terms = terms(object@fixed)
          namTerms = attr(terms, "term.labels")
          if (attr(terms, "intercept") > 0) {
              namTerms = c("(Intercept)", namTerms)
          }
          names(val$terms) = namTerms
          namTerms = factor(assign, labels = namTerms)
          attr(val, "assign") = split(order(assign), namTerms)
          val
      })
#          N <- nrow(X)
#          p <- ncol(X)
#          Qp1 <- Q + 1
#          namX <- colnames(X)
#          ngrps <- rev(ngrps)[-(1:2)]
#          stratNam <- c(names(ngrps), "Residual")
#          dfX <- dfTerms <- c(ngrps, N) - c(0, ngrps)
#          names(dfX) <- names(dfTerms) <- stratNam
#          valX <- double(p)
#          names(valX) <- namX
#          namTerms <- names(assign)
#          valTerms <- double(length(assign))
#          names(valTerms) <- namTerms
#          if (any(notIntX <- apply(X, 2, function(el) any(el != el[1])))) {
#              ## percentage of groups for which columns of X are inner
#              innP <- .Call("nlme_inner_perc_table",
#                            object,
#                            PACKAGE = "lme4")
#              dimnames(innP) <- list(namX, stratNam)
#              ## strata in which columns of X are estimated
#              ## ignoring fractional inner percentages for now
#              stratX <- stratNam[apply(innP, 1, function(el, index) max(index[el > 0]),
#                                       index = 1:Qp1)]
#              ## strata in which terms are estimated
#              notIntTerms <- unlist(lapply(assign,
#                                           function(el, notIntX) {
#                                               any(notIntX[el])
#                                           }, notIntX = notIntX))
#              stratTerms <- stratNam[unlist(lapply(assign,
#                                                   function(el, stratX, stratNam) {
#                                                       max(match(stratX[el], stratNam))
#                                                   },
#                                                   stratX = stratX, stratNam = stratNam))][notIntTerms]
#              stratX <- stratX[notIntX]
#              xDF <- table(stratX)
#              dfX[names(xDF)] <- dfX[names(xDF)] - xDF
#              if (!all(notIntX)) {                # correcting df for intercept
#                  dfX[1] <- dfX[1] - 1
#              } else {
#                  dfX[-1] <- dfX[-1] + 1
#              }
#              valX[notIntX] <- dfX[stratX]
#              ## number of parameters in each term
#              pTerms <- unlist(lapply(assign, length))[notIntTerms]
#              tDF <- tapply(pTerms, stratTerms, sum)
#              dfTerms[names(tDF)] <- dfTerms[names(tDF)] - tDF
#              if (!all(notIntTerms)) {
#                  dfTerms[1] <- dfTerms[1] - 1
#              } else {
#                  dfTerms[-1] <- dfTerms[-1] + 1
#              }
#              valTerms[notIntTerms] <- dfTerms[stratTerms]
#          } else {
#              notIntTerms <- unlist(lapply(assign,
#                                           function(el, notIntX) {
#                                               any(notIntX[el])
#                                           }, notIntX = notIntX))
#          }
#          if (!all(notIntX)) {  #intercept included
#              valX[!notIntX] <- max(dfX)
#              if (!all(notIntTerms)) {
#                  valTerms[!notIntTerms] <- max(dfTerms)
#              }
#          }
#          val <- list(X = valX, terms = valTerms)
#          attr(val, "assign") <- assign
#          val
#      })

setMethod("show", signature(object="summary.reStruct"),
          function(object) {
              digits = max(3, getOption("digits") - 2)
              useScale = length(object@useScale) > 0 && object@useScale[1]
              sc = ifelse(useScale, object@scale,  1.)
              reStdDev = lapply(object@reSumry, function(x, sc) sc*x@cor@stdDev,
                                sc = sc)
              reLens = unlist(lapply(reStdDev, length))
              reMat = array('', c(sum(reLens), 4),
                             list(rep('', sum(reLens)),
                                  c("Groups", "Name", "Variance", "Std.Dev.")))
              reMat[1+cumsum(reLens)-reLens, 1] = names(reLens)
              reMat[,2] = unlist(lapply(reStdDev, names))
              reMat[,3] = format(unlist(reStdDev)^2, digits = digits)
              reMat[,4] = format(unlist(reStdDev), digits = digits)
              if (any(reLens > 1) &&
                  !all(sapply(object@reSumry,
                              function(x) x@noCorrelation))) {
                  maxlen = max(reLens)
                  corr =
                      do.call("rbind",
                              lapply(object@reSumry,
                                     function(x, maxlen) {
                                         if (x@noCorrelation) {
                                             matrix("", dim(x@cor)[1], maxlen)
                                         } else {
                                             cc = format(round(x@cor, 3),
                                                         nsmall = 3)
                                             cc[!lower.tri(cc)] = ""
                                             nr = dim(cc)[1]
                                             cbind(cc, matrix("",
                                                              nr, maxlen-nr))
                                         }
                                     }, maxlen))
                  colnames(corr) = c("Corr", rep("", maxlen - 1))
                  reMat = cbind(reMat, corr)
              }
              if (useScale) {
                  reMat = rbind(reMat, c("Residual", "",
                                          format(sc^2, digits = digits),
                                          format(sc, digits = digits),
                                          rep('', ncol(reMat) - 4)))
              }
              cat("Random effects:\n")
              print(reMat, quote = FALSE)
              cm = object@coefficients
              if (useScale) {
                  cm[,"Std. Error"] = sc * cm[, "Std. Error"]
                  stat = cm[,1]/cm[,2]
                  pval = 2*pt(abs(stat), cm[,3], lower = FALSE)
                  nms = colnames(cm)
                  cm = cbind(cm, stat, pval)
                  colnames(cm) = c(nms, "t value", "Pr(>|t|)")
              } else {
                  cat("\nEstimated scale (compare to 1) ", object@scale, "\n")
                  stat = cm[,1]/cm[,2]
                  pval = 2*pnorm(abs(stat), lower = FALSE)
                  nms = colnames(cm)
                  cm = cbind(cm, stat, pval)
                  colnames(cm) = c(nms, "z value", "Pr(>|z|)")
              }
              cat("\nFixed effects:",
                  paste(deparse(object@fixed),
                        sep = '\n', collapse = '\n'), "\n")
              print.coefmat(cm, tst.ind = 4, zap.ind = 3)
              if (length(object@showCorrelation) > 0 && object@showCorrelation[1]) {
                  correl = object@corFixed
                  rn = rownames(cm)
                  dimnames(correl) = list(
                          abbreviate(rn, minlen=11), abbreviate(rn, minlen=6))
                  if (!is.null(correl)) {
                      p = NCOL(correl)
                      if (p > 1) {
                          cat("\nCorrelation of Fixed Effects:\n")
                          correl = format(round(correl, 3), nsmall = 3)
                          correl[!lower.tri(correl)] = ""
                          print(correl[-1, -p, drop=FALSE], quote = FALSE)
                      }
                  }
              }
          })


### Local variables:
### mode: R
### End:

setReplaceMethod("LMEoptimize", signature(x="ssclme", value="list"),
                 function(x, value)
             {
                 if (value$msMaxIter < 1) return(x)
                 st = coef(x, unconst = TRUE) # starting values
                 if (value$optimizer == "optim") {
                     optimRes =
                         if (value$analyticGradient) {
                             optim(st,
                                   fn = function(pars) {
                                       coef(x, unconst = TRUE) = pars
                                       deviance(x, REML = value$REML)
                                   },
                                   gr = function(pars) {
                                       coef(x, unconst = TRUE) = pars
                                       gradient(x, REML = value$REML,
                                                unconst = TRUE)
                                   },
                                   method = "BFGS",
                                   control = list(trace = value$msVerbose,
                                                  reltol = value$msTol,
                                                  maxit = value$msMaxIter))
                         } else {
                             optim(st,
                                   fn = function(pars) {
                                       coef(x, unconst = TRUE) = pars
                                       deviance(x, REML = value$REML)
                                   },
                                   method = "BFGS",
                                   control = list(trace = value$msVerbose,
                                                  reltol = value$msTol,
                                                  maxit = value$msMaxIter))
                         }
                     if (optimRes$convergence != 0) {
                         warning("optim failed to converge")
                     }
                     coef(x, unconst = TRUE) = optimRes$par
                 } else {
                     typsize <- rep(1.0, length(st))
                     if (is.null(value$nlmStepMax))
                         value$nlmStepMax <-
                             max(100 * sqrt(sum((st/typsize)^2)), 100)
                     nlmRes =
                         nlm(f = if (value$analyticGradient) {
                             function(pars) {
                                 coef(x, unconst = TRUE) = pars
                                 ans = deviance(x, REML = value$REML)
                                 attr(ans, "gradient") =
                                     gradient(x, REML = value$REML,
                                              unconst = TRUE)
                                 ans
                             }
                         } else {
                             function(pars)
                             {
                                 coef(x, unconst = TRUE) = pars
                                 deviance(x, REML = value$REML)
                             }
                         },
                             p = st,
                             print.level = if (value$msVerbose) 2 else 0,
                             steptol = value$msTol,
                             gradtol = value$msTol,
                             stepmax = value$nlmStepMax,
                             typsize=typsize,
                             iterlim = value$msMaxIter)
                     coef(x, unconst = TRUE) = nlmRes$estimate
                 }
                 return(x)
             })

setMethod("deviance", signature(object = "ssclme"),
          function(object, REML = FALSE, ...) {
              .Call("ssclme_factor", object, PACKAGE = "Matrix")
              object@deviance[ifelse(REML, 2, 1)]
          })

setMethod("ranef", signature(object = "ssclme"),
          function(object, ...) {
              val = .Call("ssclme_ranef", object, PACKAGE = "Matrix")
              bv = object@bVar
              names(val) = names(bv)
              for (i in seq(along = val)) {
                  dimnames(val[[i]]) = dimnames(bv[[i]])[-1]
              }
              lapply(val, t)
          })

setMethod("fixef", signature(object = "ssclme"),
          function(object, ...) {
              val = .Call("ssclme_fixef", object, PACKAGE = "Matrix")
              names(val) = dimnames(object@XtX)[[2]][seq(along = val)]
              val
          })

setMethod("vcov", signature(object = "ssclme"),
          function(object, REML = TRUE, useScale = TRUE,...) {
              ## force an "ssclme_invert"
              sc = .Call("ssclme_sigma", object, REML, PACKAGE = "Matrix")
              rr = object@RXX
              nr = nrow(rr)
              rr = rr[-nr, -nr, drop = FALSE]
              rr = rr %*% t(rr)
              if (useScale) {
                  rr = sc^2 * rr
              }
              rr
          })

setMethod("VarCorr", signature(x = "ssclme"),
          function(x, REML = TRUE, useScale = TRUE, ...) {
              val = .Call("ssclme_variances", x, PACKAGE = "Matrix")
              bVar = x@bVar
              for (i in seq(along = val)) {
                  dimnames(val[[i]]) = dimnames(bVar[[i]])[1:2]
                  val[[i]] = as(as(val[[i]], "pdmatrix"), "corrmatrix")
              }
              new("VarCorr",
                  scale = .Call("ssclme_sigma", x, REML),
                  reSumry = val,
                  useScale = useScale)
          })

setMethod("gradient", signature(x = "ssclme"),
          function(x, REML, unconst, ...)
          .Call("ssclme_gradient", x, REML, unconst))

setMethod("summary", "ssclme",
          function(object, REML = TRUE, useScale = TRUE, ...) {
              fcoef = fixef(object)
              corF <- as(as(vcov(object, REML, useScale), "pdmatrix"),
                         "corrmatrix")
              DF = getFixDF(object)
              coefs = cbind(fcoef, corF@stdDev, DF)
              nc = object@nc
              dimnames(coefs) =
                  list(names(fcoef), c("Estimate", "Std. Error", "DF"))
              ngrps = diff(object@Gp)
              ngrps = as.integer(ngrps/nc[seq(a = ngrps)])
              names(ngrps) = names(object@Omega)
              new("summary.ssclme",
                  coefficients = as.matrix(coefs),
                  scale = .Call("ssclme_sigma", object, REML),
                  denomDF = as.integer(DF),
                  REML = REML,
                  ngrps = ngrps,
                  nobs = nc[length(nc)],
                  corFixed = corF,
                  VarCorr = VarCorr(object, REML, useScale),
                  useScale = useScale,
                  showCorrelation = FALSE)
          })

## calculates degrees of freedom for fixed effects Wald tests
## This is a placeholder.  The answers are generally wrong.  It will
## be very tricky to decide what a 'right' answer should be with
## crossed random effects.

setMethod("getFixDF", signature(object="ssclme"),
          function(object, ...)
      {
          nc = object@nc[-seq(along = object@Omega)]
          p = nc[1] - 1
          n = nc[2]
          rep(n-p, p)
      })

#          Q <- length(object@random)-2
#          columns = object@random[["*fixed*"]]@columns
#          X = object@original[, columns, drop = FALSE]
#          ngrps = unlist(lapply(object@random, function(lmeLevel)
#                                 lmeLevel@nlev))
#          names(ngrps) = names(object@random)
#          val = .Call("nlme_getFixDF", object, PACKAGE = "lme4")
#          names(val$X) =
#              colnames(object@original)[object@random[["*fixed*"]]@columns]
#                                        # Convert R's assign to S-PLUS style
#          assign = object@assign.X
#          terms = terms(object@fixed)
#          namTerms = attr(terms, "term.labels")
#          if (attr(terms, "intercept") > 0) {
#              namTerms = c("(Intercept)", namTerms)
#          }
#          names(val$terms) = namTerms
#          namTerms = factor(assign, labels = namTerms)
#          attr(val, "assign") = split(order(assign), namTerms)
#          val
#       })

setMethod("show", signature(object="summary.ssclme"),
          function(object) {
              digits = max(3, getOption("digits") - 2)
              useScale = length(object@useScale) > 0 && object@useScale[1]
              REML = length(object@REML) > 0 && object@REML[1]
              cat("Random effects:\n")
              show(object@VarCorr)
              cat(sprintf("# of obs: %d, groups: ", object@nobs))
              cat(paste(paste(names(object@ngrps), object@ngrps, sep = ", "), collapse = "; "))
              cat("\n")
              if (!useScale)
                  cat("\nEstimated scale (compare to 1) ", object@scale, "\n")
              cm = object@coefficients
              if (nrow(cm) > 0) {
                  if (useScale) {
                      stat = cm[,1]/cm[,2]
                      pval = 2*pt(abs(stat), cm[,3], lower = FALSE)
                      nms = colnames(cm)
                      cm = cbind(cm, stat, pval)
                      colnames(cm) = c(nms, "t value", "Pr(>|t|)")
                  } else {
                      cm = cm[, 1:2, drop = FALSE]
                      stat = cm[,1]/cm[,2]
                      pval = 2*pnorm(abs(stat), lower = FALSE)
                      nms = colnames(cm)
                      cm = cbind(cm, stat, pval)
                      colnames(cm) = c(nms, "z value", "Pr(>|z|)")
                  }
                  cat("\nFixed effects:\n")
                  printCoefmat(cm, tst.ind = 4, zap.ind = 3)
                  if (length(object@showCorrelation) > 0 &&
                      object@showCorrelation[1]) {
                      correl = object@corFixed
                      rn = rownames(cm)
                      dimnames(correl) = list(
                              abbreviate(rn, minlen=11),
                              abbreviate(rn, minlen=6))
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
              }
          })


setMethod("coef", signature(object="summary.ssclme"),
          function(object, ...)
      {
          digits = max(3, getOption("digits") - 2)
          useScale = length(object@useScale) > 0 && object@useScale[1]
          cm = object@coefficients
          if (nrow(cm) > 0) {
              if (useScale) {
                  stat = cm[,1]/cm[,2]
                  pval = 2*pt(abs(stat), cm[,3], lower = FALSE)
                  nms = colnames(cm)
                  cm = cbind(cm, stat, pval)
                  colnames(cm) = c(nms, "t value", "Pr(>|t|)")
              } else {
                  cm = cm[, 1:2, drop = FALSE]
                  stat = cm[,1]/cm[,2]
                  pval = 2*pnorm(abs(stat), lower = FALSE)
                  nms = colnames(cm)
                  cm = cbind(cm, stat, pval)
                  colnames(cm) = c(nms, "z value", "Pr(>|z|)")
              }
              printCoefmat(cm, tst.ind = 4, zap.ind = 3)
          }
      })

setMethod("anova", signature(object="ssclme"),
          function(object, ...)
      {
          foo <- object
          foo@status["factored"] <- FALSE
          .Call("ssclme_factor", foo, PACKAGE="Matrix")
          dfr <- getFixDF(foo)
          ss <- foo@RXX[ , ".response"]^2
          ssr <- ss[[".response"]]
          ss <- ss[seq(along = dfr)]
          # FIXME: This only gives single degree of freedom tests
          ms <- ss
          df <- rep(1, length(ms))
          f <- ms/(ssr/dfr)
          P <- pf(f, df, dfr, lower.tail = FALSE)
          table <- data.frame(df, ss, ms, dfr, f, P)
          dimnames(table) <-
              list(names(ss),
                   c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
          if (any(match(names(ss), "(Intercept)", nomatch = 0)))
              table <- table[-1,]
          attr(table, "heading") = "Analysis of Variance Table"
          class(table) <- c("anova", "data.frame")
          table
      })

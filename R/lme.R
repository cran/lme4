facshuffle = function(sslm, facs)       # unexported utility
{
    s2 = sslm[[2]]
    if (all(s2 == (seq(a = s2) - 1))) return(facs)
    if (getOption("verbose")) cat(" Non-trivial permutation\n")
    s1 = sslm[[1]]
    lens = diff(s1@Gp)
    lens = lens/(s1@nc[seq(a = lens)])
    ff = vector("list", length(facs))
    for (i in seq(along = lens)) {
        sq = seq(lens[i])
        perm = 1 + s2[sq]
        s2 = s2[-sq] - lens[i]
        fi = facs[[i]]
        fip = factor(perm[as.integer(fi)])
        levels(fip)[perm] = levels(fi)
        ff[[i]] = fip
    }
    ff
}

make.mf.call = function(mf, frm, random) #unexported utility
{
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    form = frm
    form[[3]] <- (~a+b)[[2]]
    form[[3]][[2]] <- frm[[3]]
    form[[3]][[3]] <-
        as.formula((parse(text=paste("~",
                          paste(names(random),
                                collapse = "+")))[[1]]))[[2]]
    for (pdm in random) {
        tmp <- form
        tmp[[3]] <- (~a+b)[[2]]
        tmp[[3]][[2]] <- form[[3]]
        tmp[[3]][[3]] <- formula(pdm)[[2]]
        form <- tmp
    }
    environment(form) = environment(formula)
    mf$formula = form
    mf$drop.unused.levels = TRUE
    mf
}

lmeControl <-                            # Control parameters for lme
  function(maxIter = 50,
           msMaxIter = 50,
           tolerance = sqrt((.Machine$double.eps)),
           niterEM = 20,
           msTol = sqrt(.Machine$double.eps),
           msScale,
           msVerbose = getOption("verbose"),
           PQLmaxIt = 20,
           .relStep = (.Machine$double.eps)^(1/3),
           nlmStepMax = NULL,
           optimizer="nlm",
           EMverbose = getOption("verbose"),
           analyticGradient = TRUE,
           analyticHessian=FALSE)
{
    if (missing(msScale)) msScale = function(start) {
        scale <- abs(start)
        nonzero <- scale > 0
        if (any(nonzero)) {
            scale[nonzero] <- 1/scale[nonzero]
            scale[!nonzero] <- median(scale[nonzero])
        }
        else {
            scale <- rep(1, length(scale))
        }
        scale
    }
    list(maxIter = maxIter,
         msMaxIter = msMaxIter,
         tolerance = tolerance,
         niterEM = niterEM,
         msTol = msTol,
         msScale = msScale,
         msVerbose = msVerbose,
         PQLmaxIt = PQLmaxIt,
         .relStep = .relStep,
         nlmStepMax = nlmStepMax,
         optimizer=optimizer,
         EMverbose=EMverbose,
         analyticHessian=analyticHessian,
         analyticGradient=analyticGradient)
}

setMethod("lme", signature(formula = "missing"),
          function(formula, data, random,
                   method = c("REML", "ML"),
                   control = list(),
                   subset, weights, na.action, offset,
                   model = TRUE, x = FALSE, y = FALSE,...)
      {
          nCall = mCall = match.call()
          resp = getResponseFormula(data)[[2]]
          cov = getCovariateFormula(data)[[2]]
          nCall$formula = eval(substitute(resp ~ cov))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "Matrix")
      })

setMethod("lme", signature(formula = "formula", data = "groupedData",
                           random = "missing"),
          function(formula, data, random,
                   method = c("REML", "ML"),
                   control = list(),
                   subset, weights, na.action, offset,
                   model = TRUE, x = FALSE, y = FALSE,...)
      {
          nCall = mCall = match.call()
          cov = formula[[3]]
          grps = getGroupsFormula(data)[[2]]
          nCall$random = eval(substitute(~ cov | grps))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "Matrix")
      })

setMethod("lme", signature(random = "formula"),
          function(formula, data, random,
                   method = c("REML", "ML"),
                   control = list(),
                   subset, weights, na.action, offset,
                   model = TRUE, x = FALSE, y = FALSE,...)
      {
          nCall = mCall = match.call()
          cov = getCovariateFormula(random)
          nms = all.vars(getGroupsFormula(random))
          lst = lapply(nms, function(f) cov)
          names(lst) = nms
          nCall$random = lst
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "Matrix")
      })

setMethod("lme", signature(formula = "formula", data = "groupedData",
                           random = "list"),
          function(formula, data, random,
                   method = c("REML", "ML"),
                   control = list(),
                   subset, weights, na.action, offset,
                   model = TRUE, x = FALSE, y = FALSE,...)
      {
          nCall = mCall = match.call()
          nCall$data <- data@data
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "Matrix")
      })
          
setMethod("lme", signature(formula = "formula",
                           random = "list"),
          function(formula, data, random,
                   method = c("REML", "ML"),
                   control = list(),
                   subset, weights, na.action, offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
          method = match.arg(method)
          random = lapply(random, formula) # formula function, not argument
          controlvals = do.call("lmeControl", control)
          controlvals$REML = method == "REML"
          data <- eval(make.mf.call(match.call(expand.dots = FALSE),
                                    formula, random), parent.frame())
          facs = lapply(names(random),
                         function(x) eval(as.name(x), envir = data))
          names(facs) = names(random)
          facs =                        # order by decreasing number of levels
              facs[rev(order(sapply(facs, function(fac)
                                    length(levels(fac)))))]
          mmats <- c(lapply(random,
                            function(x) model.matrix(formula(x), data = data)),
                     list(.Xy = cbind(model.matrix(formula, data = data),
                          .response = model.response(data))))
          obj = .Call("ssclme_create", facs, sapply(mmats, ncol),
                       PACKAGE = "Matrix")
          facs = facshuffle(obj, facs)
          obj = obj[[1]]
          .Call("ssclme_update_mm", obj, facs, mmats, PACKAGE="Matrix")
          .Call("ssclme_initial", obj, PACKAGE="Matrix")
          .Call("ssclme_EMsteps", obj, controlvals$niterEM,
                controlvals$REML, controlvals$EMverbose, PACKAGE = "Matrix")
          LMEoptimize(obj) = controlvals
          fitted = .Call("ssclme_fitted", obj, facs, mmats, TRUE, PACKAGE = "Matrix")
          residuals = mmats$.Xy[,".response"] - fitted
          if (as.logical(x)[1]) x = mmats else x = list()
          rm(mmats)
          .Call("ssclme_to_lme", match.call(), facs, x,
                if(model) data else data.frame(list()),
                method == "REML", obj, fitted, residuals, PACKAGE = "Matrix")
      })

setMethod("fitted", signature(object="lme"),
          function(object, ...)
      {
          object@fitted
      })

setMethod("residuals", signature(object="lme"),
          function(object, ...) object@residuals )

setMethod("logLik", signature(object="lme"),
          function(object, REML = object@REML, ...) {
              val = -deviance(object@rep, REML = REML)/2
              rr = object@rep
              nc = rr@nc[-seq(a = rr@Omega)]
              attr(val, "nall") = attr(val, "nobs") = nc[2]
              attr(val, "df") = nc[1] + length(coef(rr))
              attr(val, "REML") = REML 
              class(val) <- "logLik"
              val
          })

setMethod("deviance", signature(object="lme"),
          function(object, REML, ...)
          deviance(object@rep,
                   REML = ifelse(missing(REML), object@REML, REML))
          )

setMethod("summary", signature(object="lme"),
          function(object, ...) {
              llik <- logLik(object)
              resd <- residuals(object, type="pearson")
              if (length(resd) > 5) {
                  resd <- quantile(resd)
                  names(resd) <- c("Min","Q1","Med","Q3","Max")
              }
              new("summary.lme",
                  call = object@call,
                  logLik = llik,
                  re = summary(object@rep, REML = object@REML,
                               useScale = TRUE),
                  residuals = resd)
          })

setMethod("show", signature(object = "summary.lme"),
          function(object)
      {
          rdig <- 5
          cat("Linear mixed-effects model fit by ")
          cat(ifelse(object@re@REML, "REML\n", "maximum likelihood\n") )
          if (!is.null(object@call$formula)) {
              cat("Fixed:", deparse(object@call$formula),"\n")
          }
          if (!is.null(object@call$data)) {
              cat(" Data:", deparse(object@call$data), "\n")
          }
          if (!is.null(object@call$subset)) {
              cat(" Subset:",
                  deparse(asOneSidedFormula(object@call$subset)[[2]]),"\n")
          }
          llik = object@logLik
          print(data.frame(AIC = AIC(llik), BIC = BIC(llik),
                           logLik = c(object@logLik), row.names = ""))
          cat("\n")
          object@re@useScale = TRUE
          object@re@showCorrelation = TRUE
          show(object@re)
          invisible(object)
      })

setMethod("coef", signature(object = "summary.lme"),
          function(object, ...)
      {
          coef(object@re)
      })

setMethod("show", signature(object = "lme"),
          function(object)
      {
          sumry = summary(object)
          rdig <- 5
          cat("Linear mixed-effects model\n")
          if (!is.null(object@call$formula)) {
              cat("Fixed:", deparse(object@call$formula),"\n")
          }
          if (!is.null(object@call$data)) {
              cat(" Data:", deparse( object@call$data ), "\n")
          }
          if (!is.null(object@call$subset)) {
              cat(" Subset:",
                  deparse(asOneSidedFormula(object@call$subset)[[2]]),"\n")
          }
          cat(paste(" log-", ifelse(object@REML, "restricted-", ""),
                    "likelihood: ", sep = ''), logLik(object), "\n")
          sumry@re@useScale = TRUE
          sumry@re@showCorrelation = FALSE
          saveopt = options(show.signif.stars=FALSE)
          on.exit(options(saveopt))
          show(sumry@re)
          invisible(object)
      })

setMethod("anova", signature(object = "lme"),
          function(object, ...)
      {
          mCall <- match.call(expand.dots = TRUE)
          dots <- list(...)
          modp <- logical(0)
          if (length(dots))
              modp <- sapply(dots, inherits, "lme") | sapply(dots, inherits, "lm")
          if (!any(modp)) {             # only one model - use terms
              mCall$object <- substitute(object@rep)
              return(eval(mCall, parent.frame()))
          }
          opts <- dots[!modp]
          mods <- c(list(object), dots[modp])
          names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)], as.character)
          mods <- mods[order(sapply(lapply(mods, logLik, REML = FALSE), attr, "df"))]
          calls <- lapply(mods, slot, "call")
          data <- lapply(calls, "[[", "data")
          if (any(data != data[[1]])) stop("all models must be fit to the same data object")
          header <- paste("Data:", data[[1]])
          subset <- lapply(calls, "[[", "subset")
          if (any(subset != subset[[1]])) stop("all models must use the same subset")
          if (!is.null(subset[[1]]))
              header <-
                  c(header,
#                    paste("Subset", deparse(subset[[1]], forDisplay = TRUE), sep = ": "))
                    paste("Subset", deparse(subset[[1]]), sep = ": "))
          llks <- lapply(mods, logLik, REML = FALSE)
          Df <- sapply(llks, attr, "df")
          llk <- unlist(llks)
          chisq <- 2 * pmax(0, c(NA, diff(llk)))
          dfChisq <- c(NA, diff(Df))
          val <- data.frame(Df = Df,
                            AIC = sapply(llks, AIC),
                            BIC = sapply(llks, BIC),
                            logLik = llk,
                            "Chisq" = chisq,
                            "Chi Df" = dfChisq,
                            "Pr(>Chisq)" = pchisq(chisq, dfChisq, lower = FALSE),
                            check.names = FALSE)
          class(val) <- c("anova", class(val))
          attr(val, "heading") <-
              c(header, "", "Models: <fixed>: <random>",
                paste(names(mods),
#                      unlist(lapply(lapply(calls, "[[", "formula"), deparse, forDisplay = TRUE)),
#                      unlist(lapply(lapply(calls, "[[", "random"), deparse, forDisplay = TRUE)),
                      unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
                      unlist(lapply(lapply(calls, "[[", "random"), deparse)),
                      sep = ": "),"")
          val
      })

setMethod("fixef", signature(object = "lme"),
          function(object, ...)
      {
          object = object@rep
          callGeneric()
      })

setMethod("formula", "lme", function(x, ...) x@call$formula)

setMethod("plot", signature(x = "lme"),
          function(x, y, ...)
          cat("plot method for lme not yet implemented\n"))

setMethod("ranef", signature(object = "lme"),
          function(object, ...)
      {
          object = object@rep
          lapply(callGeneric(), as.data.frame)
      })

setMethod("coef", signature(object = "lme"),
          function(object, ...)
      {
          mCall <- match.call()
          mCall[[1]] <- as.name("ranef")
          random <- eval(mCall, parent.frame())
          random
      })

setMethod("update", signature(object = "lme"),
          function(object, formula., ..., evaluate = TRUE)
      {
          call <- object@call
          if (is.null(call))
              stop("need an object with call component")
          extras <- match.call(expand.dots = FALSE)$...
          if (!missing(formula.))
              call$formula <- update.formula(formula(object), formula.)
          if (length(extras) > 0) {
              existing <- !is.na(match(names(extras), names(call)))
              for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
              if (any(!existing)) {
                  call <- c(as.list(call), extras[!existing])
                  call <- as.call(call)
              }
          }
          if (evaluate)
              eval(call, parent.frame())
          else call
      })

setMethod("vcov", signature(object = "lme"),
          function(object, ...) {
              object = object@rep
              callGeneric()
          })

setMethod("VarCorr", signature(x = "lme"),
          function(x, ...) {
              x = x@rep
              callGeneric()
          })

setMethod("gradient", signature(x = "lme"),
          function(x, ...) {
              x = x@rep
              callGeneric()
          })

setMethod("confint", signature(object = "lme"),
          function (object, parm, level = 0.95, ...) 
      {
          cf <- fixef(object)
          pnames <- names(cf)
          if (missing(parm)) 
              parm <- seq(along = pnames)
          else if (is.character(parm)) 
              parm <- match(parm, pnames, nomatch = 0)
          a <- (1 - level)/2
          a <- c(a, 1 - a)
          pct <- paste(round(100 * a, 1), "%")
          ci <- array(NA, dim = c(length(parm), 2),
                      dimnames = list(pnames[parm], pct))
          ses <- sqrt(diag(vcov(object)))[parm]
          ci[] <- cf[parm] + ses * t(outer(a, getFixDF(object@rep)[parm], qt))
          ci
      })


#lme1Control <-                            # Control parameters for lme
#  function(maxIter = 50,
#           msMaxIter = 50,
#           tolerance = sqrt((.Machine$double.eps)),
#           niterEM = 20,
#           msTol = sqrt(.Machine$double.eps),
#           msScale,
#           msVerbose = getOption("verbose"),
#           PQLmaxIt = 20,
#           .relStep = (.Machine$double.eps)^(1/3),
#           nlmStepMax = NULL,
#           optimizer="nlm",
#           EMverbose = getOption("verbose"),
#           analyticGradient = TRUE,
#           analyticHessian=FALSE)
#{
#    if (missing(msScale)) msScale = function(start) {
#        scale <- abs(start)
#        nonzero <- scale > 0
#        if (any(nonzero)) {
#            scale[nonzero] <- 1/scale[nonzero]
#            scale[!nonzero] <- median(scale[nonzero])
#        }
#        else {
#            scale <- rep(1, length(scale))
#        }
#        scale
#    }
#    list(maxIter = maxIter,
#         msMaxIter = msMaxIter,
#         tolerance = tolerance,
#         niterEM = niterEM,
#         msTol = msTol,
#         msScale = msScale,
#         msVerbose = msVerbose,
#         PQLmaxIt = PQLmaxIt,
#         .relStep = .relStep,
#         nlmStepMax = nlmStepMax,
#         optimizer=optimizer,
#         EMverbose=EMverbose,
#         analyticHessian=analyticHessian,
#         analyticGradient=analyticGradient)
#}

#setMethod("lme1", signature(formula = "missing"),
#          function(formula, data, random,
#                   method = c("REML", "ML"),
#                   control = list(),
#                   subset, weights, na.action, offset,
#                   model = TRUE, x = FALSE, y = FALSE,...)
#      {
#          nCall = mCall = match.call()
#          resp = getResponseFormula(data)[[2]]
#          cov = getCovariateFormula(data)[[2]]
#          nCall$formula = eval(substitute(resp ~ cov))
#          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
#                mCall, PACKAGE = "Matrix")
#      })

#setMethod("lme1", signature(formula = "formula", data = "groupedData",
#                           random = "missing"),
#          function(formula, data, random,
#                   method = c("REML", "ML"),
#                   control = list(),
#                   subset, weights, na.action, offset,
#                   model = TRUE, x = FALSE, y = FALSE,...)
#      {
#          nCall = mCall = match.call()
#          cov = formula[[3]]
#          grps = getGroupsFormula(data)[[2]]
#          nCall$random = eval(substitute(~ cov | grps))
#          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
#                mCall, PACKAGE = "Matrix")
#      })

#setMethod("lme1", signature(random = "formula"),
#          function(formula, data, random,
#                   method = c("REML", "ML"),
#                   control = list(),
#                   subset, weights, na.action, offset,
#                   model = TRUE, x = FALSE, y = FALSE,...)
#      {
#          nCall = mCall = match.call()
#          cov = getCovariateFormula(random)
#          nms = all.vars(getGroupsFormula(random))
#          lst = lapply(nms, function(f) cov)
#          names(lst) = nms
#          nCall$random = lst
#          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
#                mCall, PACKAGE = "Matrix")
#      })

#setMethod("lme1", signature(formula = "formula", data = "groupedData",
#                           random = "list"),
#          function(formula, data, random,
#                   method = c("REML", "ML"),
#                   control = list(),
#                   subset, weights, na.action, offset,
#                   model = TRUE, x = FALSE, y = FALSE,...)
#      {
#          nCall <- mCall <- match.call()
#          nCall$data <- data@data
#          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
#                mCall, PACKAGE = "Matrix")
#      })
          
#setMethod("lme1", signature(formula = "formula",
#                           random = "list"),
#          function(formula, data, random,
#                   method = c("REML", "ML"),
#                   control = list(),
#                   subset, weights, na.action, offset,
#                   model = TRUE, x = FALSE, y = FALSE, ...)
#      {
#          method <- match.arg(method)
#          random <- lapply(random, formula) # formula function, not the argument
#          controlvals <- do.call("lme1Control", control)
#          controlvals$REML <- method == "REML"
#          data <- eval(make.mf.call(match.call(expand.dots = FALSE),
#                                           formula, random), parent.frame())
#          facs <- lapply(names(random),
#                         function(x) eval(as.name(x), envir = data))
#          names(facs) <- names(random)
#          facs <-                        # order by decreasing number of levels
#              facs[rev(order(sapply(facs, function(fac)
#                                    length(levels(fac)))))]
#          mmats <- c(lapply(random,
#                            function(x) model.matrix(formula(x), data = data)),
#                     list(.Xy = cbind(model.matrix(formula, data = data),
#                          .response = model.response(data))))
#          obj <- .Call("lmeRep_create", facs, sapply(mmats, ncol),
#                       PACKAGE = "Matrix")
#          .Call("lmeRep_update_mm", obj, facs, mmats, PACKAGE="Matrix")
#          .Call("lmeRep_initial", obj, PACKAGE="Matrix")
#          .Call("lmeRep_ECMEsteps", obj, 
#                controlvals$niterEM,
#                controlvals$REML,
#                controlvals$EMverbose,
#                TRUE,
#                PACKAGE = "Matrix")
#          LMEoptimize(obj) <- controlvals
#          #fitted = .Call("ssclme_fitted", obj, facs, mmats, TRUE, PACKAGE = "Matrix")
#          #residuals = mmats$.Xy[,".response"] - fitted
#          #if (as.logical(x)[1]) x = mmats else x = list()
#          #rm(mmats)
#          obj
#      })

#setReplaceMethod("LMEoptimize", signature(x="lmeRep", value="list"),
#                 function(x, value)
#             {
#                 if (value$msMaxIter < 1) return(x)
#                 st = coef(x, unconst = TRUE) # starting values
#                 if (value$optimizer == "optim") {
#                     optimRes =
#                         if (value$analyticGradient) {
#                             optim(st,
#                                   fn = function(pars) {
#                                       coef(x, unconst = TRUE) = pars
#                                       deviance(x, REML = value$REML)
#                                   },
#                                   gr = function(pars) {
#                                       coef(x, unconst = TRUE) = pars
#                                       gradient(x, REML = value$REML,
#                                                unconst = TRUE)
#                                   },
#                                   method = "BFGS",
#                                   control = list(trace = value$msVerbose,
#                                                  reltol = value$msTol,
#                                                  maxit = value$msMaxIter))
#                         } else {
#                             optim(st,
#                                   fn = function(pars) {
#                                       coef(x, unconst = TRUE) = pars
#                                       deviance(x, REML = value$REML)
#                                   },
#                                   method = "BFGS",
#                                   control = list(trace = value$msVerbose,
#                                                  reltol = value$msTol,
#                                                  maxit = value$msMaxIter))
#                         }
#                     if (optimRes$convergence != 0) {
#                         warning("optim failed to converge")
#                     }
#                     coef(x, unconst = TRUE) = optimRes$par
#                 } else {
#                     typsize <- rep(1.0, length(st))
#                     if (is.null(value$nlmStepMax))
#                         value$nlmStepMax <-
#                             max(100 * sqrt(sum((st/typsize)^2)), 100)
#                     nlmRes =
#                         nlm(f = if (value$analyticGradient) {
#                             function(pars) {
#                                 coef(x, unconst = TRUE) = pars
#                                 ans = deviance(x, REML = value$REML)
#                                 attr(ans, "gradient") =
#                                     gradient(x, REML = value$REML,
#                                              unconst = TRUE)
#                                 ans
#                             }
#                         } else {
#                             function(pars)
#                             {
#                                 coef(x, unconst = TRUE) = pars
#                                 deviance(x, REML = value$REML)
#                             }
#                         },
#                             p = st,
#                             print.level = if (value$msVerbose) 2 else 0,
#                             steptol = value$msTol,
#                             gradtol = value$msTol,
#                             stepmax = value$nlmStepMax,
#                             typsize=typsize,
#                             iterlim = value$msMaxIter)
#                     coef(x, unconst = TRUE) = nlmRes$estimate
#                 }
#                 return(x)
#             })

#setMethod("deviance", signature(object = "lmeRep"),
#          function(object, REML = FALSE, ...) {
#              .Call("lmeRep_factor", object, PACKAGE = "Matrix")
#              object@deviance[ifelse(REML, 2, 1)]
#          })

#setMethod("ranef", signature(object = "lmeRep"),
#          function(object, ...) {
#              .Call("lmeRep_ranef", object, PACKAGE = "Matrix")
#          })

#setMethod("fixef", signature(object = "lmeRep"),
#          function(object, ...) {
#              val = .Call("lmeRep_fixef", object, PACKAGE = "Matrix")
#              names(val) = object@cnames[[".fixed"]]
#              val[-length(val)]
#          })

#setMethod("vcov", signature(object = "lmeRep"),
#          function(object, REML = TRUE, useScale = TRUE,...) {
#              ## force an "lmeRep_invert"
#              sc <- .Call("lmeRep_sigma", object, PACKAGE = "Matrix")
#              rr <- object@RXX
#              nms <- object@cnames[[".fixed"]]
#              dimnames(rr) <- list(nms, nms)
#              nr <- nrow(rr)
#              rr <- rr[-nr, -nr, drop = FALSE]
#              rr <- rr %*% t(rr)
#              if (useScale) {
#                  rr = sc^2 * rr
#              }
#              rr
#          })

#setMethod("VarCorr", signature(x = "lmeRep"),
#          function(x, REML = TRUE, useScale = TRUE, ...) {
#              val = .Call("lmeRep_variances", x, PACKAGE = "Matrix")
#              bVar = x@bVar
#              for (i in seq(along = val)) {
#                  dimnames(val[[i]]) = dimnames(bVar[[i]])[1:2]
#                  val[[i]] = as(as(val[[i]], "pdmatrix"), "corrmatrix")
#              }
#              new("VarCorr",
#                  scale = .Call("lmeRep_sigma", x, REML),
#                  reSumry = val,
#                  useScale = useScale)
#          })

#setMethod("gradient", signature(x = "lmeRep"),
#          function(x, REML, unconst, ...)
#          .Call("lmeRep_gradient", x, REML, unconst))

#setMethod("summary", "lmeRep",
#          function(object, REML = TRUE, useScale = TRUE, ...) {
#              fcoef = fixef(object)
#              corF <- as(as(vcov(object, REML, useScale), "pdmatrix"),
#                         "corrmatrix")
#              DF = getFixDF(object)
#              coefs = cbind(fcoef, corF@stdDev, DF)
#              nc = object@nc
#              dimnames(coefs) =
#                  list(names(fcoef), c("Estimate", "Std. Error", "DF"))
#              ngrps = diff(object@Gp)
#              ngrps = as.integer(ngrps/nc[seq(a = ngrps)])
#              names(ngrps) = names(object@Omega)
#              new("summary.ssclme",
#                  coefficients = as.matrix(coefs),
#                  scale = .Call("lmeRep_sigma", object, REML),
#                  denomDF = as.integer(DF),
#                  REML = REML,
#                  ngrps = ngrps,
#                  nobs = nc[length(nc)],
#                  corFixed = corF,
#                  VarCorr = VarCorr(object, REML, useScale),
#                  useScale = useScale,
#                  showCorrelation = FALSE)
#          })

### calculates degrees of freedom for fixed effects Wald tests
### This is a placeholder.  The answers are generally wrong.  It will
### be very tricky to decide what a 'right' answer should be with
### crossed random effects.

#setMethod("getFixDF", signature(object="lmeRep"),
#          function(object, ...)
#      {
#          nc = object@nc[-seq(along = object@Omega)]
#          p = nc[1] - 1
#          n = nc[2]
#          rep(n-p, p)
#      })

#setMethod("anova", signature(object="lmeRep"),
#          function(object, ...)
#      {
#          foo <- object
#          foo@status["factored"] <- FALSE
#          .Call("lmeRep_factor", foo, PACKAGE="Matrix")
#          dfr <- getFixDF(foo)
#          ss <- foo@RXX[ , ".response"]^2
#          ssr <- ss[[".response"]]
#          ss <- ss[seq(along = dfr)]
#          names(ss) <- object@cnames[[".fixed"]][seq(along = dfr)]
#          # FIXME: This only gives single degree of freedom tests
#          ms <- ss
#          df <- rep(1, length(ms))
#          f <- ms/(ssr/dfr)
#          P <- pf(f, df, dfr, lower.tail = FALSE)
#          table <- data.frame(df, ss, ms, dfr, f, P)
#          dimnames(table) <-
#              list(names(ss),
#                   c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
#          if (any(match(names(ss), "(Intercept)", nomatch = 0)))
#              table <- table[-1,]
#          attr(table, "heading") = "Analysis of Variance Table"
#          class(table) <- c("anova", "data.frame")
#          table
#      })

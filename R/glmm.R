setMethod("GLMM", signature(data = "missing"),
          function(formula, family, data, random, ...)
      {
          nCall = mCall = match.call()
          nCall$data = list()
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call", mCall)
      })

setMethod("GLMM", signature(formula = "missing", data = "groupedData"),
          function(formula, family, data, random, ...)
      {
          nCall = mCall = match.call()
          resp = getResponseFormula(data)[[2]]
          cov = getCovariateFormula(data)[[2]]
          nCall$formula = eval(substitute(resp ~ cov))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call", mCall)
      })

setMethod("GLMM", signature(formula = "formula", data = "groupedData",
                           random = "missing"),
          function(formula, family, data, random, ...)
      {
          nCall = mCall = match.call()
          cov = formula[[3]]
          grps = getGroupsFormula(data)[[2]]
          nCall$random = eval(substitute(~ cov | grps))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call", mCall)
      })


setMethod("GLMM", signature(formula = "formula", random = "formula"),
          function(formula, family, data, random, ...)
      {
          nCall = mCall = match.call()
          nCall$random = lapply(getGroupsFormula(random, asList = TRUE),
                                function(x, form) form,
                                form = pdLogChol(getCovariateFormula(random)))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call", mCall)
      })

setMethod("GLMM", signature(formula = "formula", random = "list"),
          function(formula, family, data, random, correlation, weights,
                   control, niter, method, verbose, nEM.IRLS,
                   model, x, ...)
      {
          if (missing(nEM.IRLS))
              nEM.IRLS <- 1
          if (missing(verbose)) verbose = FALSE
          m <- Call <- match.call()
          method <- if(missing(method)) "PQL" else
                    match.arg(method, c("PQL", "Laplace"))
          nm <- names(m)[-1]
          dontkeep <-
              is.element(nm, c("correlation", "control", "niter",
                               "verbose", "nEM.IRLS", "method"))
          for(i in nm[dontkeep]) m[[i]] <- NULL
          m[[1]] <- as.name("glmmStruct")
          m$nextraCols <- 1
          m$method <- method
          fit <- eval(m, parent.frame())

          off <- fit@offset
          w <-  fit@prior.weights
          origy <- fit@origy
          fam <- fit@family

          ## We always do the PQL steps. For that, we extract the reStruct
          ## slot from fit
          fit <- as(fit, "reStruct")
          control = if (missing(control)) lmeControl() else
                    do.call("lmeControl", control)
          EMsteps(fit) <- control

          control$niterEM <- nEM.IRLS
          converged <- FALSE
          eta <- .Call("nlme_reStruct_fitted", fit, PACKAGE="lme4")

          for(i in seq(length=if(missing(niter)) 20 else niter)) {
              ##update zz and wz
              mu <- fam$linkinv(eta)
              mu.eta.val <- fam$mu.eta(eta)
              zz <- eta + (origy - mu)/mu.eta.val  - off
              wz <- w * mu.eta.val^2 / fam$variance(mu)

              response(fit) <- zz
              weighted(fit) <- sqrt(abs(wz))
              EMsteps(fit) <- control
              LMEoptimize(fit) <- control
              if(verbose) {
                  cat("iteration", i, "\n")
###             class(fit) <- "glmmStruct"
###             fit@logLik <- as.numeric(NA)
###             cat("Approximate logLik:",
###                 .Call("nlme_glmmLa2_logLikelihood",
###                     fit, NULL, PACKAGE="lme4"), "\n")
###            class(fit) <- "reStruct"
###            fit@logLik <- as.numeric(NA)
                  cat("Parameters:", coef(fit), "\n")
                  cat("Fixed Effectsf:", fixef(fit), "\n")
              }
              etaold <- eta
              eta <- .Call("nlme_reStruct_fitted", fit, PACKAGE="lme4")
              if(sum((eta-etaold)^2) < 1e-6*sum(eta^2)) {
                  converged <- TRUE
                  break
              }
          }
          if (control$msMaxIter > 0 && !converged)
              stop("IRLS iterations in GLMM failed to converge")
          ## We recreate the glmm object and set the reStruct slot
          .Call("nlme_replaceSlot", fit, "logLik",
                as.numeric(NA))
          .Call("nlme_replaceSlot", fit, "dontCopy", TRUE)
          fit <- .Call("nlme_replaceSlot", eval(m, parent.frame()),
                       "reStruct", fit)
          .Call("nlme_replaceSlot", fit, c("reStruct", "dontCopy"), TRUE)
          fit <- .Call("nlme_glmmLaplace_solveOnly", fit,
                       500, 1, PACKAGE="lme4")
          .Call("nlme_replaceSlot", fit, "reStruct",
                .Call("nlme_commonDecompose", fit@reStruct, NULL, PACKAGE="lme4"))
          .Call("nlme_replaceSlot", fit, c("reStruct", "dontCopy"), FALSE)

          if (method != "PQL") {
              ## Do the 2nd order Laplace fit here
          }
          ## zero some of the matrix slots
          if (!missing(x) && x == FALSE)
              .Call("nlme_replaceSlot", fit, c("reStruct", "original"),
                    matrix(0.0, nrow = 0, ncol = 0))
          .Call("nlme_replaceSlot", fit, c("reStruct", "decomposed"),
                matrix(0.0, nrow = 0, ncol = 0))
          .Call("nlme_replaceSlot", fit, c("reStruct", "weighted"),
                matrix(0.0, nrow = 0, ncol = 0))

          if (!missing(model) && model == FALSE)
              .Call("nlme_replaceSlot", fit, "frame",
                    data.frame())
          .Call("nlme_replaceSlot", fit, "fitted",
                fam$linkinv(if (is.null(fit@na.action)) {
                    fitted(fit@reStruct)[fit@reStruct@reverseOrder]
                } else {
                    napredict(attr(data, "na.action"),
                              fitted(fit@reStruct)[fit@reStruct@reverseOrder])
                }))
          fit
      })

setMethod("summary", signature(object="glmm"),
          function(object, ...) {
              llik <- logLik(object)    # has an oldClass
              resd <- residuals(object, type="pearson")
              if (length(resd) > 5) {
                  resd <- quantile(resd)
                  names(resd) <- c("Min","Q1","Med","Q3","Max")
              }
              ans <- new("summary.glmm",
                         call = object@call,
                         logLik = llik,
                         AIC = AIC(llik),
                         BIC = BIC(llik),
                         re = summary(as(object, "reStruct")),
                         residuals = resd,
                         method = object@method,
                         family = object@family$family,
                         link = object@family$link)
              ans@re@useScale = !(ans@family %in% c("binomial", "poisson"))
              ans
          })

setMethod("show", "summary.glmm",
          function(object) {
              rdig <- 5
              cat("Generalized linear mixed-effects model fit by ")
              cat(switch(object@method, PQL="PQL\n",
                         Laplace="2nd order Laplace\n"))
              cat(" Family:", object@family, "with",
                  object@link, "link\n")
              cat(" Data:", deparse( object@call$data ), "\n")
              if (!is.null(object@call$subset)) {
                  cat("  Subset:",
                      deparse(asOneSidedFormula(object@call$subset)[[2]]),"\n")
              }
              print(data.frame(AIC = object@AIC, BIC = object@BIC,
                               logLik = c(object@logLik), row.names = ""))
              cat("\n")
              object@re@showCorrelation = TRUE
              show(object@re)
              ## Should this be part of the show method for summary.reStruct?
              cat("\nNumber of Observations:", object@re@nobs)
              cat("\nNumber of Groups: ")
              ngrps <- object@re@ngrps
              if ((length(ngrps)) == 1) {
                  cat(ngrps,"\n")
              } else {				# multiple nesting
                  cat("\n")
                  print(ngrps)
              }
              invisible(object)
          })

setMethod("show", "glmm",
          function(object)
      {
          sumry = summary(object)
          rdig <- 5
          cat("Generalized linear mixed-effects model\n")
          cat(" Family:", sumry@family, "with",
              sumry@link, "link\n")
          cat(" Data:", deparse( sumry@call$data ), "\n")
          if (!is.null(sumry@call$subset)) {
              cat("  Subset:",
                  deparse(asOneSidedFormula(sumry@call$subset)[[2]]),"\n")
          }
          cat(paste(" log-", ifelse(sumry@re@REML, "restricted-", ""),
                    "likelihood: ", sep = ''), sumry@logLik, "\n")
          sumry@re@showCorrelation = FALSE
          saveopt = options(show.signif.stars=FALSE)
          on.exit(saveopt)
          show(sumry@re)
          options(saveopt)
          on.exit()
          cat("\nNumber of Observations:", sumry@re@nobs, "\n")
          invisible(object)
      })

setMethod("getResponse", signature(object="glmm"),
          function(object, form)
      {
          object@origy
      })

setMethod("logLik", signature(object="glmm"),
          function(object)
      {
          value = .Call("nlme_glmmLaplace_logLikelihood", object,
                         NULL,          # do not pass new parameter value
                         50, 1, PACKAGE="lme4")
          p = length(object@reStruct@random[["*fixed*"]]@columns)
          # df calculated from sigma + fixed effects + random effects pars
          attr(value, "df") = 1 + p +
              sum(unlist(lapply(object@reStruct@random,
                                function(x)length(coef(x)))))
          attr(value, "nall") = length(object@fitted)
          attr(value, "nobs") = length(object@fitted) - p
          class(value) = "logLik"
          value
      })

### Local variables:
### mode: R
### End:

lmeControl <-
  ## Control parameters for lme
  function(maxIter = 50, msMaxIter = 50, tolerance =
           sqrt((.Machine$double.eps)), niterEM = 35,
           msTol = sqrt(.Machine$double.eps), msScale, msVerbose = FALSE,
           returnObject = FALSE, gradHess = TRUE, apVar = TRUE,
           .relStep = (.Machine$double.eps)^(1/3), minAbsParApVar = 0.05,
           nlmStepMax = NULL,
           natural = TRUE, optimizer="nlm", EMverbose=FALSE,
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
    list(maxIter = maxIter, msMaxIter = msMaxIter, tolerance = tolerance,
         niterEM = niterEM, msTol = msTol, msScale = msScale,
         msVerbose = msVerbose, returnObject = returnObject,
         gradHess = gradHess , apVar = apVar, .relStep = .relStep,
         nlmStepMax = nlmStepMax,
         minAbsParApVar = minAbsParApVar, natural = natural,
         optimizer=optimizer, EMverbose=EMverbose,
         analyticHessian=analyticHessian,
         analyticGradient=analyticGradient)
}

setMethod("lme", signature(data = "missing"),
          function(formula, data, random, correlation, weights, subset,
                   method, na.action, control, model, x)
      {
          nCall = mCall = match.call()
          nCall$data = list()
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "lme4")
      })

setMethod("lme", signature(formula = "missing", data = "groupedData"),
          function(formula, data, random, correlation, weights, subset,
                   method, na.action, control, model, x)
      {
          nCall = mCall = match.call()
          resp = getResponseFormula(data)[[2]]
          cov = getCovariateFormula(data)[[2]]
          nCall$formula = eval(substitute(resp ~ cov))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "lme4")
      })

setMethod("lme", signature(formula = "formula", data = "groupedData",
                           random = "missing"),
          function(formula, data, random, correlation, weights, subset,
                   method, na.action, control, model, x)
      {
          nCall = mCall = match.call()
          cov = formula[[3]]
          grps = getGroupsFormula(data)[[2]]
          nCall$random = eval(substitute(~ cov | grps))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "lme4")
      })


setMethod("lme", signature(formula = "formula", random = "formula"),
          function(formula, data, random, correlation, weights, subset,
                   method, na.action, control, model, x)
      {
          nCall = mCall = match.call()
          nCall$random = lapply(getGroupsFormula(random, asList = TRUE),
                                function(x, form) form,
                                form = pdLogChol(getCovariateFormula(random)))

          nCall$data <- as(data, "data.frame")

          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "lme4")
      })

setMethod("lme", signature(formula = "formula", random = "list"),
          function(formula, data, random, correlation, weights, subset,
                   method, na.action, control, model, x)
      {
          if (missing(model))
              model = TRUE
          if (missing(x))
              x = TRUE
          random = lapply(random, function(x)
                          if(inherits(x, "formula")) pdLogChol(x) else x)
          method = if (missing(method)) "REML" else
                   match.arg(method, c("REML", "ML"))
          controlvals <- if (missing(control)) lmeControl() else
                            do.call("lmeControl", control)
          mCall <- match.call(expand.dots = FALSE)
          mCall[[1]] <- as.name("model.frame")
          names(mCall)[2] <- "formula"
          mCall$random <- mCall$correlation <- mCall$method <-
              mCall$control <- NULL
          form <- formula
          form[[3]] <- (~a+b)[[2]]
          form[[3]][[2]] <- formula[[3]]
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
          environment(form) <- environment(formula)
          mCall$formula <- form
          mCall$drop.unused.levels <- TRUE
          data <- eval(mCall, parent.frame())
          re <- reStruct(fixed = formula, random = random,
                         data = data,
                         REML = method != "ML",
                         analyticHessian=controlvals$analyticHessian)
          .Call("nlme_replaceSlot", re, "dontCopy", TRUE, PACKAGE = "lme4")
          .Call("nlme_replaceSlot", re, "analyticHessian",
                FALSE, PACKAGE = "lme4")
          EMsteps(re) <- controlvals
          .Call("nlme_replaceSlot", re, "analyticHessian",
                controlvals$analyticHessian, PACKAGE = "lme4")
          LMEoptimize(re) <- controlvals
          .Call("nlme_replaceSlot", re, "dontCopy", FALSE, PACKAGE = "lme4")
          ## zero some of the matrix slots
          if (x == FALSE)
              .Call("nlme_replaceSlot", re, "original",
                    matrix(0.0, nrow = 0, ncol = 0), PACKAGE = "lme4")
          .Call("nlme_replaceSlot", re, "decomposed",
                matrix(0.0, nrow = 0, ncol = 0), PACKAGE = "lme4")
          .Call("nlme_replaceSlot", re, "weighted",
                matrix(0.0, nrow = 0, ncol = 0), PACKAGE = "lme4")
          if (model == FALSE)
              data = data.frame()
          new("lme", reStruct = re, call = match.call(),
              fitted = if (is.null(attr(data, "na.action"))) {
                  fitted(re)[re@reverseOrder]
              } else {
                  napredict(attr(data, "na.action"), fitted(re)[re@reverseOrder])
              },
              frame = data, na.action = attr(data, "na.action"))
      })

setAs("lme", "reStruct",
      function(from) from@reStruct,
      function(from, value) {
          if (nrow(from@frame) != nrow(value@original))
              stop("Dimension mismatch between model.frame and original matrix")
          from@reStruct <- value
          from
      })

setMethod("fitted", signature=c(object="lme"),
          function(object, ...)
      {
          object@fitted
      })


setMethod("residuals", signature=c(object="lme"),
          function(object, ...)
      {
          re <- as(object, "reStruct")
          if (is.null(object@na.action)) {
              (getResponse(object@reStruct) - fitted(object@reStruct))[re@reverseOrder]
          } else {
              napredict(object@na.action,
                        (getResponse(object@reStruct) -
                         fitted(object@reStruct))[object@reStruct@reverseOrder])
          }
      })


setMethod("logLik", signature(object="lme"),
          function(object) logLik(object@reStruct))

setMethod("summary", signature(object="lme"),
          function(object, ...) {
              llik <- logLik(object)    # has an oldClass
              resd <- residuals(object, type="pearson")
              if (length(resd) > 5) {
                  resd <- quantile(resd)
                  names(resd) <- c("Min","Q1","Med","Q3","Max")
              }
              new("summary.lme",
                  call = object@call,
                  logLik = llik,
                  AIC = AIC(llik),
                  BIC = BIC(llik),
                  re = summary(as(object, "reStruct")),
                  residuals = resd)
          })

setMethod("show", "summary.lme",
          function(object) {
              rdig <- 5
              cat("Linear mixed-effects model fit by ")
              cat(ifelse(object@re@REML, "REML\n", "maximum likelihood\n") )
              cat(" Data:", deparse( object@call$data ), "\n")
              if (!is.null(object@call$subset)) {
                  cat("  Subset:",
                      deparse(asOneSidedFormula(object@call$subset)[[2]]),"\n")
              }
              print(data.frame(AIC = object@AIC, BIC = object@BIC,
                               logLik = c(object@logLik), row.names = ""))
              cat("\n")
              object@re@useScale = TRUE
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

setMethod("show", "lme",
          function(object)
      {
          sumry = summary(object)
          rdig <- 5
          cat("Linear mixed-effects model\n")
          cat(" Data:", deparse( sumry@call$data ), "\n")
          if (!is.null(sumry@call$subset)) {
              cat("  Subset:",
                  deparse(asOneSidedFormula(sumry@call$subset)[[2]]),"\n")
          }
          cat(paste(" log-", ifelse(sumry@re@REML, "restricted-", ""),
                    "likelihood: ", sep = ''), sumry@logLik, "\n")
          sumry@re@useScale = TRUE
          sumry@re@showCorrelation = FALSE
          saveopt = options(show.signif.stars=FALSE)
          on.exit(saveopt)
          show(sumry@re)
          options(saveopt)
          on.exit()
          cat("\nNumber of Observations:", sumry@re@nobs, "\n")
          invisible(object)
      })


setMethod("isInitialized", "lmeLevelList",
          function(object) all(sapply(object[seq(length=length(object)-2)],
                                      function(x) isInitialized(x@precision))),
          valueClass = "logical")

setMethod("anova", signature(object = "lme"),
          function(object, ...)
          cat("anova method for lme not yet implemented\n"))

setMethod("fixef", signature(object = "lme"),
          function(object, ...)
      {
          object = object@reStruct
          callGeneric()
      })

setMethod("formula", "lme", function(x, ...) x@call$formula)

setMethod("intervals", signature(object = "lme", level = "ANY"),
          function(object, level = 0.95, ...)
          cat("intervals method for lme not yet implemented\n"))

setMethod("plot", signature(x = "lme"),
          function(x, y, ...)
          cat("plot method for lme not yet implemented\n"))

setMethod("ranef", signature(object = "lme"),
          function(object, ...)
      {
          object = object@reStruct
          callGeneric()
      })

setMethod("coef", signature(object = "lme"),
          function(object, ...)
      {
          object = object@reStruct
          callGeneric()
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

setMethod("getGroups", signature(object="reStruct",
                                 form="missing",
                                 data="missing",
                                 sep="missing"),
          function(object, form, level, data, sep)
      {
          object <- object@reStruct
          callGeneric()
      })

setMethod("getResponse", signature(object="lme"),
          function(object, form)
      {
          object <- object@reStruct
          callGeneric()
      })

### Local variables:
### mode: R
### End:

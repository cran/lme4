

          

setMethod("lmList", signature(formula = "formula", data = "data.frame"),
          function(formula, data, level, subset, na.action, pool)
      {
          mCall = frmCall = match.call()
          resp = getResponseFormula(formula)[[2]]
          cov = getCovariateFormula(formula)[[2]]
          lmForm = eval(substitute(resp ~ cov))
          gfm = getGroupsFormula(formula)
          if (is.null(gfm)) gfm = getGroupsFormula(data)
          if (is.null(gfm))
              stop("Unable to determine a grouping formula from either the formula or the data")
          val <- lapply(split(data, eval(gfm[[2]], data)),
                        function(dat, formula)
                    {
                        ans <- try({
                            data <- as.data.frame(dat)
                            lm(formula = formula, data = data)
                        })
                        if (inherits(ans, "try-error"))
                            NULL
                        else ans
                    }, formula = lmForm)
          if (missing(pool)) pool = TRUE
          new("lmList", val, call = mCall, pool = pool)
      })


setMethod("coef", signature(object = "lmList"),
            ## Extract the coefficients and form a  data.frame if possible
          function(object, augFrame = FALSE, data = NULL,
                   which = NULL, FUN = mean, omitGroupingFactor = TRUE, ...)
      {
          coefs = lapply(object, coef)
          non.null = !unlist(lapply(coefs, is.null))
          if (sum(non.null) > 0) {
              template = coefs[non.null][[1]]
              if (is.numeric(template)) {
                  co <- matrix(template,
                               ncol = length(template),
                               nrow = length(coefs),
                               byrow = TRUE,
                               dimnames = list(names(object), names(template)))
                  for (i in names(object)) {
                      co[i,] = if (is.null(coefs[[i]])) { NA } else coefs[[i]]
                  }
                  coefs = as.data.frame(co)
                  effectNames = names(coefs)
                  if(augFrame) {
                      if (is.null(data)) {
                          data = getData(object)
                      }
                      data = as.data.frame(data)
                      if (is.null(which)) {
                          which = 1:ncol(data)
                      }
                      data = data[, which, drop = FALSE]
                      ## eliminating columns with same names as effects
                      data = data[, is.na(match(names(data), effectNames)), drop = FALSE]
                      data = gsummary(data, FUN = FUN, groups = getGroups(object))
                      if (omitGroupingFactor) {
                          data <- data[, is.na(match(names(data),
                                                     names(getGroupsFormula(object,
                                                                            asList = TRUE)))),
                                       drop = FALSE]
                      }
                      if (length(data) > 0) {
                          coefs = cbind(coefs, data[row.names(coefs),,drop = FALSE])
                      }
                  }
                  attr(coefs, "level") = attr(object, "level")
                  attr(coefs, "label") = "Coefficients"
                  attr(coefs, "effectNames") = effectNames
                  attr(coefs, "standardized") = FALSE
                  #attr(coefs, "grpNames") <- deparse(getGroupsFormula(object)[[2]])
                  #class(coefs) <- c("coef.lmList", "ranef.lmList", class(coefs))
              }
          }
          coefs
      })

setMethod("show", signature(object = "lmList"),
          function(object)
      {
          mCall = object@call
          cat("Call:", deparse(mCall), "\n")
          cat("Coefficients:\n")
          invisible(print(coef(object)))
          if (object@pool) {
              cat("\n")
              poolSD = pooledSD(object)
              dfRes = attr(poolSD, "df")
              RSE = c(poolSD)
              cat("Degrees of freedom: ", length(unlist(lapply(object, fitted))),
                  " total; ", dfRes, " residual\n", sep = "")
              cat("Residual standard error:", format(RSE))
              cat("\n")
          }
      })

setMethod("pooledSD", signature(object = "lmList"),
          function(object)
      {
          sumsqr <- apply(sapply(object,
                                 function(el) {
                                     if (is.null(el)) {
                                         c(0,0)
                                     } else {
                                         res = resid(el)
                                         c(sum(res^2), length(res) - length(coef(el)))
                                  }
                              }), 1, sum)
          if (sumsqr[2] == 0) {
              stop("No degrees of freedom for estimating std. dev.")
          }
          val <- sqrt(sumsqr[1]/sumsqr[2])
          attr(val, "df") <- sumsqr[2]
          val
      })

setMethod("intervals", signature(object = "lmList", level = "ANY"),
          function(object, level = 0.95, ...)
          cat("intervals method for lmList not yet implemented\n"))

setMethod("plot", signature(x = "lmList"),
          function(x, y, ...)
          cat("plot method for lmList not yet implemented\n"))

setMethod("update", signature(object = "lmList"),
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

setMethod("formula", signature(x = "lmList"),
          function(x, ...) x@call[["formula"]])

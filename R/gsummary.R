setMethod("gsummary", signature(groups = "missing"),
          function (object, FUN, form, level, groups,
                    omitGroupingFactor = FALSE, 
                    invariantsOnly = FALSE, ...)
      {
          gCall <- nCall <- match.call()
          gCall[[1]] <- as.name("getGroups")
          nCall$groups <- eval(gCall, parent.frame())
          eval(nCall, parent.frame())
      })

setMethod("gsummary",
          signature(object = "groupedData", groups = "factor"),
          function (object, FUN, form, level, groups,
                    omitGroupingFactor = FALSE, 
                    invariantsOnly = FALSE, ...)
      {
          nCall <- match.call()
          nCall$object <- substitute(object@data, list(object = nCall$object))
          eval(nCall, parent.frame())
      })
          
setMethod("gsummary", signature(object = "data.frame", groups = "factor"),
          function (object, FUN, form, level, groups,
                    omitGroupingFactor = FALSE, 
                    invariantsOnly = FALSE, ...)
      {
          if (missing(FUN)) FUN <- function(x) mean(x, na.rm = TRUE)
          groups <- groups[drop = TRUE]
          gunique <- levels(groups)
          firstInGroup <- match(gunique, groups)
          asFirst <- firstInGroup[match(groups, gunique)]
          value <- as.data.frame(object[firstInGroup, , drop = FALSE])
          row.names(value) <- as.character(gunique)
          varying <-
              unlist(lapply(object,
                            function(x, first) any(x != x[first]),
                            first = asFirst))
          if (any(varying) && (!invariantsOnly)) {
              Mode <- function(x) names(which.max(table(x)))
              if (is(FUN, "function")) {
                  FUN <- list(numeric = FUN, ordered = Mode, factor = Mode)
              }
              else {
                  if (!(is.list(FUN) &&
                        all(sapply(FUN, class) == "function")))
                      stop("FUN can only be a function or a list of functions")
                  FF <- FUN
                  FUN <- list(numeric = function(x) mean(x, na.rm = TRUE),
                              ordered = Mode, factor = Mode)
                  FUN[names(FF)] <- FF
              }
              for (nm in names(object)[varying]) {
                  dClass <- data.class(object[[nm]])
                  if (dClass == "numeric") {
                      value[[nm]] <- as.vector(tapply(object[[nm]], 
                                                      groups, FUN[["numeric"]], ...))
                  }
                  else {
                      value[[nm]] <- as.vector(tapply(as.character(object[[nm]]), 
                                                      groups, FUN[[dClass]]))
                      if (inherits(object[, nm], "ordered")) {
                          value[[nm]] <- ordered(value[, nm], levels = levels(object[, 
                                                              nm]))[drop = TRUE]
                      }
                      else {
                          value[[nm]] <- factor(value[, nm], levels = levels(object[, 
                                                             nm]))[drop = TRUE]
                      }
                  }
              }
          }
          else {
              value <- value[, !varying, drop = FALSE]
          }
          if (omitGroupingFactor) {
              if (is.null(form)) {
                  stop("Cannot omit grouping factor without \"form\"")
              }
              grpForm <- getGroupsFormula(form, asList = TRUE)
              if (missing(level)) 
                  level <- length(grpForm)
              grpNames <- names(grpForm)[level]
              whichKeep <- is.na(match(names(value), grpNames))
              if (any(whichKeep)) {
                  value <- value[, whichKeep, drop = FALSE]
              }
              else {
                  return(NULL)
              }
          }
          value
      })

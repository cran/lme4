# Provisional method.
# The definition of the method should change when groupedData is defined
#  as an S4 class. (update: which has now happened)

# setMethod("formula",
#           "groupedData",
#           function(x, ...) attr(x, "formula"),
#           valueClass = "formula")



## experimental groupedData class -- deepayan


groupedData <-
    function(formula, data = sys.parent(1),  ## parent.frame() ? is this a data.frame ?
             order.groups = TRUE,
             FUN = function(x) max(x, na.rm = TRUE),
             outer = ~0, inner = ~0,
             labels = list(), units = list())
{
    new("groupedData",
        data = data,
        formula = formula,
        outer = outer,
        inner = inner,
        labels = labels,
        units = units)
}



setMethod("show", "groupedData",
          function(object) {
              cat("groupedData with Formula: ")
              show(formula(object))
              cat("\n")
              show(object@data)
          })
          


convertFromS3groupedData <- function(from)
{
    data <- from
    class(data) <- "data.frame"
    attr(data, "labels") <- NULL
    attr(data, "units") <- NULL
    attr(data, "inner") <- NULL
    attr(data, "outer") <- NULL
    attr(data, "formula") <- NULL
    attr(data, "order.groups") <- NULL

    ans <- groupedData(data = data,
                       formula = attr(from, "formula"),
                       outer = if (is.null(attr(from, "outer")))
                       formula(~0) else attr(from, "outer"),
                       inner = if (is.null(attr(from, "inner")))
                       formula(~0) else attr(from, "inner"),
                       labels = if (is.null(attr(from, "labels")))
                       list() else attr(from, "labels"),
                       units = if (is.null(attr(from, "units")))
                       list() else attr(from, "units"))
    ans
}

## Method for x[i, j]

setMethod("[",
          signature(x = "groupedData", i = "missing", j = "ANY",
                    drop = "logical"),

          function(x, i, j, drop) {
              if (drop) warning("drop=TRUE will be ignored")
              x@data <- x@data[, j, drop = FALSE]
              x
          })


setMethod("[",
          signature(x = "groupedData", i = "ANY", j = "missing",
                    drop = "logical"),

          function(x, i, j, drop) {
              if (drop) warning("drop=TRUE will be ignored")
              x@data <-
                  if (nargs() == 4) x@data[i, , drop = FALSE]
                  else x@data[i, drop = FALSE]
              x
          })


setMethod("[",
          signature(x = "groupedData", i = "missing", j = "ANY",
                    drop = "missing"),

          function(x, i, j, drop) {
              x@data <- x@data[, j, drop = FALSE]
              x
          })


setMethod("[",
          signature(x = "groupedData", i = "ANY", j = "missing",
                    drop = "missing"),

          function(x, i, j, drop) {
              x@data <-
                  if (nargs() == 3) x@data[i, ]
                  else x@data[i]
              x
          })


setMethod("[",
          signature(x = "groupedData", i = "ANY", j = "ANY",
                    drop = "missing"),

          function(x, i, j, drop) {
              x@data <- x@data[i, j, drop = FALSE]
              x
          })


setMethod("[",
          signature(x = "groupedData", i = "ANY", j = "ANY",
                    drop = "logical"),
          
          function(x, i, j, drop) {
              if (drop) warning("drop=TRUE will be ignored")
              x@data <- x@data[i, j, drop = FALSE]
              x
          })


## Method for x[[i, j]]


setMethod("[[",
          signature(x = "groupedData", i = "ANY", j = "missing"),
          function(x, i, j) {
              x@data[[i]]
          })

setMethod("[[",
          signature(x = "groupedData", i = "ANY", j = "ANY"),
          function(x, i, j) {
              x@data[[i, j]]
          })


## Method for x$name

setMethod("$",
          signature(x = "groupedData", name = "ANY"),
          function(x, name) {
              x@data[[name]]
          })



## Method for x$name <-


setMethod("$<-",
          signature(x = "groupedData", name = "ANY", value = "ANY"),
          function(x, name, value) {
              x@data[[name]] <- value
              x
          })





## Method for x[...] <-


setMethod("[<-",
          signature(x = "groupedData", i = "ANY", j = "ANY", value = "ANY"),
          function(x, i, j, value) {
              x@data[i, j] <- value
              x
          })

setMethod("[<-",
          signature(x = "groupedData", i = "ANY", j = "missing", value = "ANY"),
          function(x, i, j, value) {
              if (nargs() == 3) x@data[i] <- value
              else x@data[i,] <- value
              x
          })

setMethod("[<-",
          signature(x = "groupedData", i = "missing", j = "ANY", value = "ANY"),
          function(x, i, j, value) {
              x@data[, j] <- value
              x
          })
















setMethod("formula",
          "groupedData",
          function(x, ...) x@formula)

setAs("groupedData", "data.frame",
      function(from) from@data)

collectorPlot <-
    function(formula, data, displayfunction,
             labels, xargs,
             ...)
{
    ## It's probably sufficient to have displayfunction a quoted fn name (like "xyplot")
    ## Think about it
    require("lattice", quietly = TRUE)
    dotargs <- list(...)


    if (!is.function(displayfunction)) {
        ## needs more sophistication/checking ?
        displayfunction <- get(displayfunction)
    }

    ## adjust arguments - dotargs and xargs

    # add labels
    xargs$xlab <- labels$x
    xargs$ylab <- labels$y

    ## Combine args:
    dotargs <- c(dotargs, xargs[!(names(xargs) %in% names(dotargs))])

    do.call("displayfunction",
            c(list(formula = formula, data = data), dotargs))
}

setMethod("summary", signature(object = "groupedData"),
          function(object, ...) summary(as(object, "data.frame"), ...))

setMethod("plot",
          signature(x = "groupedData", y="missing"),
          function(x,
                   formula,
                   inner = NULL, outer = NULL,
                   ...)
      {
          respVar <- getResponseFormula(x)[[2]]
          isNumericResponse <- is.numeric(eval(respVar, x@data))
          
          covVar <- getCovariateFormula(x)[[2]]
          if (!(isTrivialCovariate <- covVar == 1))
              isNumericCovariate <- is.numeric(eval(covVar, x@data))
          
          groupVar <- getGroupsFormula(x)[[2]]
          isMultipleGroups <- length(groupVar) > 2
          if (isMultipleGroups)
              return("multiple groups not implemented yet")
          if (!isNumericResponse)
              return("non-numeric response plots not implemented yet")

          if (isTrivialCovariate || !isNumericCovariate) {
              xargs <- list()
              displayfunction = "dotplot"
              form = as.formula(paste(groupVar, respVar, sep = "~"))
              groups = if (isTrivialCovariate) NULL else covVar
              
              xlab <-
                  if ("y" %in% names(x@labels)) x@labels$y
                  else as.character(respVar)
              if ("y" %in% names(x@units)) xlab <- paste(xlab, x@units$y)
              ylab <- as.character(groupVar)
              labels <- list(x = xlab, y = ylab)

              if (!is.null(groups))
                  xargs$groups <- groups

              ## groups needs to be combined with inner, perhaps
              ##if (!isTrivialCovariate && is.null(inner)) inner = covVar
          } else if (isNumericCovariate) {
              ## Custom panel function
              panel.nfn = function(x, y, ...) {
                  panel.xyplot(x, y, ...)
                  y.avg <- tapply(y, x, mean) # lines through average y
                  y.avg <- y.avg[!is.na(y.avg)]
                  if (length(y.avg) > 0) {
                      xvals <- as.numeric(names(y.avg))
                      ord <- order(xvals)
                      panel.xyplot(xvals[ord], y.avg[ord], type = "l", ...)
                  }
              }

              xlab <-
                  if ("x" %in% names(x@labels)) x@labels$x
                  else as.character(covVar)
              if ("x" %in% names(x@units)) xlab <- paste(xlab, x@units$x)
              ylab <-
                  if ("y" %in% names(x@labels)) x@labels$y
                  else as.character(respVar)
              if ("y" %in% names(x@units)) ylab <- paste(ylab, x@units$y)
              labels <- list(x = xlab, y = ylab)

              xargs <- list(aspect = "xy", grid = TRUE)

              displayfunction = "xyplot"
              form = as.formula( paste(   paste(respVar, covVar, sep = "~"), groupVar, sep = "|") )
              
              groups = NULL
              if (is.null(inner)) inner <- x@inner
              if (inherits(inner, "formula") && !inner[[2]]==0) groups <- inner[[2]]
              ## what's innerGroups ?

              if (!is.null(groups)) {
                  xargs$groups <- groups
                  xargs$panel <- function(x, y, grid = TRUE, ...)
                  {
                      if (grid) panel.grid()
                      panel.superpose(x, y, ...)
                  }
                  xargs$panel.groups <- panel.nfn
              } else {
                  xargs$panel <- function(x, y, grid = TRUE, ...)
                  {
                      if (grid) panel.grid()
                      panel.nfn(x, y, ...)
                  }
              }
              
              
          }
          else return("didn't expect this to happen!")

          collectorPlot(formula = if (missing(formula)) form else formula,
                        data = x@data,
                        displayfunction = displayfunction,
                        labels = labels,
                        xargs = xargs,
                        ...)
      })

## Tools for plot methods

## was old nlme plot method (???)
plot.nmGroupedData <-
      function(x, collapseLevel = Q, displayLevel = collapseLevel,
                  outer = NULL, inner = NULL, preserve = NULL, FUN = mean,
                          subset = NULL, key = TRUE, grid = TRUE, ...)
{
      args <- list(outer = outer, inner = inner, key = key, grid = grid, ...)
        Q <- length(getGroupsFormula(x, asList = TRUE))
        if (is.null(preserve) && (collapseLevel < Q) && (!is.null(inner))) {
                if (is.logical(inner)) {
                          preserve <- attr(x, "inner")[[displayLevel]]
                      } else {
                                preserve <- inner
                            }
            }
        x <- collapse(x, collapseLevel, displayLevel, outer, inner,
                      preserve, FUN, subset)
        args[["innerGroups"]] <- attr(x, "innerGroups")
        args[["x"]] <- x
        do.call("plot", args)
  }








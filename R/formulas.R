## Methods for dealing with formulas

splitFormula <-
    function(form, sep = "/")
{
    ## split, on the sep call, the rhs of a formula into a list of subformulas
    if (inherits(form, "formula") ||
        mode(form) == "call" && form[[1]] == as.name("~"))
        return(splitFormula(form[[length(form)]], sep = sep))
    if (mode(form) == "call" && form[[1]] == as.name(sep))
        return(do.call("c", lapply(as.list(form[-1]), splitFormula, sep = sep)))
    if (mode(form) == "(") return(splitFormula(form[[2]], sep = sep))
    if (length(form) < 1) return(NULL)
    list(asOneSidedFormula(form))
}

subFormula <- function(form, pos = 2)
{
    ## extract component pos of form as a formula preserving the environment
    comp = form[[pos]]
    val = eval(substitute(~comp))
    environment(val) = environment(form)
    val
}

getCovariateFormula <- function(object)
{
    ## Return the primary covariate formula as a one sided formula
    form <- formula(object)
    form <- form[[length(form)]]
    if (length(form) == 3 && form[[1]] == as.name("|")){ # conditional expression
        form <- form[[2]]
    }
    eval(substitute(~form))
}

getResponseFormula <- function(object)
{
    ## Return the response formula as a one sided formula
    form <- formula(object)
    if (!(inherits(form, "formula") && (length(form) == 3)))
        stop("object must yield a two-sided formula")
    subFormula(form, 2)
}

setMethod("getGroupsFormula", signature(object = "ANY"),
          function(object, asList = FALSE, sep = "/")
      {
          form = formula(object)
          if (!inherits(form, "formula")) stop("object must yield a formula")
          rhs = form[[length(form)]]
          if (length(rhs) < 2 || rhs[[1]] != as.name("|")) return(NULL)
          if (asList) {
              val = splitFormula(asOneSidedFormula(rhs[[3]]), sep = sep)
              names(val) = unlist(lapply(val, function(el) deparse(el[[2]])))
              return(val)
          }
          asOneSidedFormula(rhs[[3]])
      })


setMethod("getGroups", signature(object="data.frame", form="formula"),
          function(object, form, level, data, sep)
              eval(getGroupsFormula(form)[[2]], object))

setMethod("getGroups", signature(object="groupedData"),
          function(object, form, level, data, sep) {
              object = object@data
              callGeneric()
          })


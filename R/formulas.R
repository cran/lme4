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
          function(object, form, level, data, sep, ...)
              eval(getGroupsFormula(form)[[2]], object))

# Return the pairs of expressions separated by vertical bars

findbars <- function(term)
{
    if (is.name(term) || is.numeric(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    if (!is.call(term)) stop("term must be of class call")
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}

# Return the formula omitting the pairs of expressions separated by vertical bars

nobars <- function(term)
{
    # FIXME: is the is.name in the condition redundant?
    #   A name won't satisfy the first condition.
    if (!('|' %in% all.names(term)) || is.name(term)) return(term)
    if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
    if (length(term) == 2) {
        nb <- nobars(term[[2]])
        if (is.null(nb)) return(NULL)
        term[[2]] <- nb
        return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

# Substitute the '+' function for the '|' function

subbars <- function(term)
{
    if (is.name(term) || is.numeric(term)) return(term)
    if (length(term) == 2) {
        term[[2]] <- subbars(term[[2]])
        return(term)
    }
    stopifnot(length(term) == 3)
    if (is.call(term) && term[[1]] == as.name('|')) term[[1]] <- as.name('+')
    term[[2]] <- subbars(term[[2]])
    term[[3]] <- subbars(term[[3]])
    term
}
    

              

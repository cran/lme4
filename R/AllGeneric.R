if (!isGeneric("getGroups")) {
    ## Return the groups associated with object according to form.
    setGeneric("getGroups",
               function(object, form, level, data, sep, ...)
               standardGeneric("getGroups"))
}

if (!isGeneric("getGroupsFormula")) {
    ## Return the formula(s) for the groups associated with object.
    ## The result is a one-sided formula unless asList is TRUE in which case
    ## it is a list of formulas, one for each level.
    setGeneric("getGroupsFormula",
               function(object, asList = FALSE, sep = "/")
               standardGeneric("getGroupsFormula"))
}

if (!isGeneric("getCovariate")) {
    ## Return the primary covariate associated with object
    setGeneric("getCovariate",
               function(object, form = formula(object), data = list())
               standardGeneric("getCovariate"))
}

if (!isGeneric("getResponse")) {
    ## Return the response variable associated with object
    setGeneric("getResponse",
               function(object, form = formula(object))
               standardGeneric("getResponse"))
}

setGeneric("lmer",
           function(formula, data, family,
                    method = c("REML", "ML", "PQL", "Laplace", "AGQ"),
                    control = list(),
                    subset, weights, na.action, offset,
                    model = TRUE, x = FALSE, y = FALSE,...)
           standardGeneric("lmer"))

if (!isGeneric("LMEoptimize<-")) {
    setGeneric("LMEoptimize<-", function(x, ..., value)
               standardGeneric("LMEoptimize<-"))
}

if (!isGeneric("fixef")) {
    setGeneric("fixef", function(object, ...) standardGeneric("fixef"))
}

fixed.effects <- function(object, ...) {
    ## fixed.effects was an alternative name for fixef
    .Deprecated("fixef")
    mCall = match.call()
    mCall[[1]] = as.name("fixef")
    eval(mCall, parent.frame())
}

if (!isGeneric("ranef")) {
    setGeneric("ranef", function(object, ...)
               standardGeneric("ranef"))
}

random.effects <- function(object, ...) {
    ## random.effects was an alternative name for ranef
    .Deprecated("ranef")
    mCall = match.call()
    mCall[[1]] = as.name("ranef")
    eval(mCall, parent.frame())
}

if (!isGeneric("BIC")) {
    setGeneric("BIC", function(object, ...) standardGeneric("BIC"))
}

setMethod("BIC", "logLik",
          function(object, ...)
          -2 * (c(object) - attr(object, "df") * log(attr(object, "nobs"))/2)
          )


if (!isGeneric("lmList")) {
    setGeneric("lmList",
               function(formula, data, level, subset, na.action, pool)
               standardGeneric("lmList"))
}

if (!isGeneric("pooledSD")) {
    setGeneric("pooledSD", function(object) standardGeneric("pooledSD"))
}

if (!isGeneric("VarCorr")) {
    setGeneric("VarCorr", function(x, ...) standardGeneric("VarCorr"))
}

if (!isGeneric("gradient")) {           # not exported
    setGeneric("gradient", function(x, ...) standardGeneric("gradient"))
}

if (!isGeneric("getFixDF")) {           # not exported
    setGeneric("getFixDF", function(object, ...) standardGeneric("getFixDF"))
}

if (!isGeneric("gsummary")) {
    setGeneric("gsummary",
          function (object, FUN, form, level, groups,
                    omitGroupingFactor = FALSE, 
                    invariantsOnly = FALSE, ...)
               standardGeneric("gsummary"))
}

if (!isGeneric("param")) {              # not exported
    setGeneric("param", function(object, ...) standardGeneric("param"))
}

if (!isGeneric("param<-")) {              # not exported
    setGeneric("param<-", function(object, ..., value)
               standardGeneric("param<-"))
}


if (!isGeneric("lmList")) {
    setGeneric("lmList",
               function(formula, data, family, subset, weights,
                        na.action, offset, pool, ...)
               standardGeneric("lmList"))
}

if (!isGeneric("pooledSD")) {
    setGeneric("pooledSD", function(object) standardGeneric("pooledSD"))
}

if (!isGeneric("gsummary")) {
    setGeneric("gsummary",
          function (object, FUN, form, level, groups,
                    omitGroupingFactor = FALSE, 
                    invariantsOnly = FALSE, ...)
               standardGeneric("gsummary"))
}

## ----------------------- lmer-related Generics ---------------------------

if (!isGeneric("isNested"))
    setGeneric("isNested", function(x, ...) standardGeneric("isNested"))

if (!isGeneric("LMEoptimize<-")) {
    setGeneric("LMEoptimize<-", function(x, ..., value)
               standardGeneric("LMEoptimize<-"))
}

if (!isGeneric("fixef")) {
    setGeneric("fixef", function(object, ...) standardGeneric("fixef"))
}

if (!isGeneric("denomDF")) {
    setGeneric("denomDF", function(x, ...) standardGeneric("denomDF"))
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

if (!isGeneric("VarCorr")) {
    setGeneric("VarCorr", function(x, ...) standardGeneric("VarCorr"))
}

if (!isGeneric("postVar")) {            # posterior variances
    setGeneric("postVar", function(object, ...)
               standardGeneric("postVar"))
}

if (!isGeneric("gradient")) {           # not exported
    setGeneric("gradient", function(x, ...) standardGeneric("gradient"))
}

if (!isGeneric("getFixDF")) {           # not exported
    setGeneric("getFixDF", function(object, ...) standardGeneric("getFixDF"))
}

if (!isGeneric("mcmcsamp")) {
    setGeneric("mcmcsamp", function(object, n = 1, verbose = FALSE, ...)
	       standardGeneric("mcmcsamp"))
}

if (!exists("simulate", mode = "function")) {
    setGeneric("simulate", function(object, nsim = 1, seed = NULL, ...)
               standardGeneric("simulate"))
}


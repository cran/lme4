setGeneric("lmList",
           function(formula, data, family, subset, weights,
                    na.action, offset, pool, ...)
           standardGeneric("lmList"))

setGeneric("pooledSD", function(object) standardGeneric("pooledSD"))

setGeneric("gsummary",
           function (object, FUN, form, level, groups,
                     omitGroupingFactor = FALSE, 
                     invariantsOnly = FALSE, ...)
           standardGeneric("gsummary"))

setGeneric("isNested", function(x, ...) standardGeneric("isNested"))

setGeneric("LMEoptimize<-", function(x, ..., value)
           standardGeneric("LMEoptimize<-"))

setGeneric("fixef", function(object, ...) standardGeneric("fixef"))

setGeneric("denomDF", function(x, ...) standardGeneric("denomDF"))

fixed.effects <- function(object, ...) {
    ## fixed.effects was an alternative name for fixef
    .Deprecated("fixef")
    mCall = match.call()
    mCall[[1]] = as.name("fixef")
    eval(mCall, parent.frame())
}

setGeneric("ranef", function(object, ...) standardGeneric("ranef"))

random.effects <- function(object, ...) {
    ## random.effects was an alternative name for ranef
    .Deprecated("ranef")
    mCall = match.call()
    mCall[[1]] = as.name("ranef")
    eval(mCall, parent.frame())
}

setGeneric("BIC", function(object, ...) standardGeneric("BIC"))

setMethod("BIC", "logLik",
          function(object, ...)
          -2 * (c(object) - attr(object, "df") * log(attr(object, "nobs"))/2)
          )

setGeneric("VarCorr", function(x, ...) standardGeneric("VarCorr"))

setGeneric("postVar", function(object, ...) standardGeneric("postVar"))

setGeneric("gradient", function(x, ...) standardGeneric("gradient"))

setGeneric("getFixDF", function(object, ...) standardGeneric("getFixDF"))

setGeneric("mcmcsamp", function(object, n = 1, verbose = FALSE, ...)
           standardGeneric("mcmcsamp"))

setGeneric("simulate", function(object, nsim = 1, seed = NULL, ...)
           standardGeneric("simulate"))


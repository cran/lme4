## The generics Names, and Names<- will be deprecated in nlme_4.0
if (!isGeneric("Names")) {
    setGeneric("Names", function(object, ...)
           {
               .Deprecated("names")
               standardGeneric("Names")
           })
}

if (!isGeneric("Names<-")) {
    setGeneric("Names<-",
               function(object, value)
           {
               .Deprecated("names<-")
               standardGeneric("Names<-")
           })
}

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

setGeneric("lme",
          function(formula, data, random,
                   method = c("REML", "ML"),
                   control = list(),
                   subset, weights, na.action, offset,
                   model = TRUE, x = FALSE, y = FALSE,...)
           standardGeneric("lme"))

if (!isGeneric("LMEoptimize<-")) {
    setGeneric("LMEoptimize<-", function(x, ..., value)
               standardGeneric("LMEoptimize<-"))
}

if (!isGeneric("fixef")) {
    setGeneric("fixef", function(object, ...) standardGeneric("fixef"))
}

fixed.effects = function(object, ...) {
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

random.effects = function(object, ...) {
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

## FIXME: Can this be replaced by confint?
if (!isGeneric("intervals")) {
    setGeneric("intervals",
               function(object, level = 0.95, ...)
               standardGeneric("intervals"))
}

if (!isGeneric("lmList")) {
    setGeneric("lmList",
               function(formula, data, level, subset, na.action, pool)
               standardGeneric("lmList"))
}

if (!isGeneric("GLMM")) {
    setGeneric("GLMM",
               function(formula, family, data, random,
                        method = c("PQL", "Laplace"),
                        control = list(),
                        subset,
                        weights,
                        na.action,
                        offset,
                        model = TRUE, x = FALSE, y = FALSE, ...)
               standardGeneric("GLMM"))
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

## The generics Names, Names<-, pdFactor, pdMatrix, and corMatrix
## will be deprecated in nlme_4.0
if (!isGeneric("Names")) {
    setGeneric("Names",
               function(object, ...)
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

if (!isGeneric("pdFactor")) {
    setGeneric("pdFactor",
               function(object) standardGeneric("pdFactor"))
}

if (!isGeneric("pdMatrix")) {
    setGeneric("pdMatrix",
               function(object) standardGeneric("pdMatrix"))
}

if (!isGeneric("corFactor")) {
    ## extractor for transpose inverse square root factor of corr matrix
    setGeneric("corFactor",
               function(object, ...) standardGeneric("corFactor"))
}

if (!isGeneric("corMatrix")) {
    ## extractor for correlation matrix or the transpose inverse
    ## square root matrix
    setGeneric("corMatrix",
               function(object, ...) standardGeneric("corMatrix"))
}

if (!isGeneric("isInitialized")) {
    setGeneric("isInitialized",
               function(object)
               #standardGeneric("isInitialized"), valueClass = "logical")
               standardGeneric("isInitialized"))
}

if (!isGeneric("logDet")) {
    setGeneric("logDet",
               function(object, covariate = getCovariate(object), ...)
               standardGeneric("logDet"))#,
               #valueClass = "numeric")
}

if (!isGeneric("matrix<-")) {
    setGeneric("matrix<-",
               function(object, value)
               standardGeneric("matrix<-"))
}

if (!isGeneric("coef<-")) {
    setGeneric("coef<-",
               function(object, value)
               standardGeneric("coef<-"))
}

if (!isGeneric("weighted<-")) {
    setGeneric("weighted<-", function(x, value) standardGeneric("weighted<-"))
}

if (!isGeneric("model.matrix<-")) {
    setGeneric("model.matrix<-", function(x, value)
               standardGeneric("model.matrix<-"))
}

if (!isGeneric("getGroups")) {
    ## Return the groups associated with object according to form.
    setGeneric("getGroups",
               function(object, form, level, data, sep)
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
    ## Return the primary covariate associated with object
    setGeneric("getResponse",
               function(object, form = formula(object))
               standardGeneric("getResponse"))
}

if (!isGeneric("LMEgradient")) {
    setGeneric("LMEgradient",
               function(x, A, nlev) standardGeneric("LMEgradient"))
}

setGeneric("lme",
           function(formula, data, random, correlation, weights, subset,
                    method, na.action, control, model, x)
           standardGeneric("lme"))

if (!isGeneric("EMupdate<-")) {
    setGeneric("EMupdate<-",
               function(x, nlev, value) standardGeneric("EMupdate<-"))
}

if (!isGeneric("reStruct")) {
    setGeneric("reStruct",
               function(fixed, random, data, weights, REML, nextraCols=0)
               standardGeneric("reStruct"))
}

if (!isGeneric("glmmStruct")) {
    setGeneric("glmmStruct",
               function(formula, random, family, data, nextraCols=0,
                        method="PQL", ...)
               standardGeneric("glmmStruct"))
}

if (!isGeneric("pdgradient")) {
    setGeneric("pdgradient", function(x) standardGeneric("pdgradient"))
}

if (!isGeneric("EMsteps<-")) {
    setGeneric("EMsteps<-", function(x, value) standardGeneric("EMsteps<-"))
}

if (!isGeneric("LMEoptimize<-")) {
    setGeneric("LMEoptimize<-", function(x, value)
               standardGeneric("LMEoptimize<-"))
}

if (!isGeneric("response<-")) {
    setGeneric("response<-", function(x, value)
               standardGeneric("response<-"))
}

if (!isGeneric("fixef")) {
    setGeneric("fixef", function(object, ...)
               standardGeneric("fixef"))
}

if (!isGeneric("fixef<-")) {
    setGeneric("fixef<-",
               function(x, value) standardGeneric("fixef<-"))
}

## fixed.effects was an alternative name
fixed.effects = function(object, ...) {
    .Deprecated("fixef")
    mCall = match.call()
    mCall[[1]] = as.name("fixef")
    eval(mCall, parent.frame())
}

if (!isGeneric("ranef")) {
    setGeneric("ranef", function(object, ...)
               standardGeneric("ranef"))
}

## random.effects was an alternative name for ranef
random.effects = function(object, ...) {
    .Deprecated("ranef")
    mCall = match.call()
    mCall[[1]] = as.name("ranef")
    eval(mCall, parent.frame())
}

if (!isGeneric("BIC")) {
    setGeneric("BIC", function(object, ...)
               standardGeneric("BIC"))
}

setMethod("BIC", "logLik",
          function(object, ...)
          -2 * (c(object) - attr(object, "df") * log(attr(object, "nobs"))/2)
          )

if (!isGeneric("getFixDF")) {
    ## Return the formula(s) for the groups associated with object.
    ## The result is a one-sided formula unless asList is TRUE in which case
    ## it is a list of formulas, one for each level.
    setGeneric("getFixDF", function(object) standardGeneric("getFixDF"))
}

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
               function(formula, family, data, random, control, niter,
                        method, verbose, ...)
               standardGeneric("GLMM"))
}

if (!isGeneric("pooledSD")) {
    setGeneric("pooledSD", function(object) standardGeneric("pooledSD"))
}

if (!isGeneric("VarCorr")) {
    setGeneric("VarCorr", function(x)
               standardGeneric("VarCorr"))
}

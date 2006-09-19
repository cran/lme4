## Class definitions for the package

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

setClass("lmList.confint", contains = "array")

## -------------------- lmer-related Classes --------------------------------

setOldClass("data.frame")
setOldClass("family")
setOldClass("logLik")
setOldClass("terms")

## mixed effects representation
setClass("mer",
	 representation(## original data
			flist = "list", # list of grouping factors
			Zt = "dgCMatrix",  # sparse representation of Z'
			X = "matrix",	   # X
			y = "numeric",	   # y
			wts = "numeric",   # weights
			wrkres = "numeric",# working residuals (copy of y for LMMs)
			method = "character", # parameter estimation method
			useScale = "logical", # should scale factor be included
                        family = "family",
			call = "call",	   # call to model-fitting function
			## invariants derived from data structure
			cnames = "list",   # column names of model matrices
			nc = "integer",	   # dimensions of blocks in Omega
			Gp = "integer",	   # Pointers to groups of rows in Zt
			## quantities that vary when Z, X or y are updated
			XtX = "dpoMatrix", # X'X
			ZtZ = "dsCMatrix", # Z'Z
			ZtX = "dgeMatrix", # Z'X
			Zty = "numeric",   # Z'y
			Xty = "numeric",   # X'y
			## primary slots that vary during the optimization
			## When Omega is updated, these are updated
			Omega = "list", # list of relative precision matrices
			## Cholesky factor of inflated [Z:X:y]'[Z:X:y]
			L = "dCHMsuper", # sparse Cholesky factor of Z'Z + Omega
			RZX = "dgeMatrix",
			RXX = "dtrMatrix",
			rZy = "numeric",
			rXy = "numeric",
			devComp = "numeric", # Components of deviance
			deviance = "numeric", # Current deviance (ML and REML)
			## Secondary slots only evaluated when requested.
			fixef = "numeric",
			ranef = "numeric",
			RZXinv = "dgeMatrix",
			bVar = "list",
			gradComp = "list",
			## status indicator
			status = "logical"
			)
	)

## Representation of linear and generalized linear mixed effects model
setClass("lmer",
	 representation(frame = "data.frame",
			terms = "terms"),
	 contains = "mer")

setClass("glmer",
	 representation(#family = "family", # glm family - move here later
                        frame = "data.frame",
			terms = "terms",
                        weights = "numeric"),
	 contains = "mer")

setClass("summary.mer", # the "mer" result ``enhanced'' :
	 representation(
			isG   = "logical",
			methTitle = "character",
			logLik= "logLik",
			ngrps = "integer",
			sigma = "numeric", # scale, non-negative number
			coefs = "matrix",
			vcov = "dpoMatrix",
			REmat = "matrix",
			AICtab= "data.frame"
			),
	 contains = "mer")

setClass("summary.lmer", contains = c("summary.mer", "lmer"))

setClass("ranef.lmer", contains = "list")

setClass("coef.lmer", contains = "list")

setClass("pedigree", representation =
	 list(sire = "integer", dam = "integer", label = "character"),
	 validity = function(object) {
	     n <- length(sire <- object@sire)
	     if (length(dam <- object@dam) != n)
		 return("sire and dam slots must be the same length")
	     if (length(object@label) != n)
		 return("'label' slot must have the same length as 'sire' and 'dam'")
	     if(n == 0) return(TRUE)
	     animal <- 1:n
	     snmiss <- !is.na(sire)
	     dnmiss <- !is.na(dam)
	     if (any(sire[snmiss] >= animal[snmiss]) ||
		 any(dam[dnmiss] >= animal[dnmiss]))
		 return("the sire and dam must precede the offspring")
             if (any(sire[snmiss] < 1 | sire[snmiss] > n) |
                 any(dam[dnmiss] < 1 | dam[dnmiss] > n))
                 return(paste("Non-missing sire or dam must be in [1,",
                              n, "]", sep = ''))
	     TRUE
	 })

## ----------------------- lmer-related Generics ---------------------------

## Hmm: If this does not match *exactly* the "formula" - method in ./lmer.R
## ---  the  match.call() in there may give a very different result
setGeneric("lmer",
	   function(formula, data, family = gaussian,
		    method = c("REML", "ML", "PQL", "Laplace", "AGQ"),
		    control = list(), start, subset, weights, na.action,
		    offset, contrasts = NULL, model = TRUE,
		    ...)
	   standardGeneric("lmer"))

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


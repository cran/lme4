.onLoad <- function(lib, pkg)
{
    if ("package:nlme" %in% search()) {
        stop(paste("Package lme4 conflicts with package nlme.\n",
                   "To attach lme4 you must restart R without package nlme."))
    }
}

setOldClass("data.frame")
setOldClass("family")
setOldClass("logLik")
setOldClass("terms")

setClass("lme", representation(call = "call",
                               facs = "list",
                               x = "list",
                               model = "data.frame",
                               REML = "logical",
                               rep = "ssclme",
                               fitted = "numeric",
                               residuals = "numeric",
                               terms = "terms",
                               assign = "integer"))

setClass("GLMM", representation(family = "family",
                                logLik = "numeric",
                                fixef = "numeric",
                                Hessian = "matrix",
                                method = "character"),
         contains = "lme")

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

setClass("VarCorr",
         representation(scale="numeric",
                        reSumry="list",
                        useScale="logical"),
         prototype = list(scale = 1.0, useScale = TRUE))

setClass("summary.ssclme",
         representation(coefficients="matrix",
                        scale="numeric",
                        denomDF="integer",
                        REML="logical",
                        ngrps="integer",
                        nobs="integer",
                        corFixed="corrmatrix",
                        VarCorr="VarCorr",
                        useScale="logical",
                        showCorrelation="logical"
                        ))

setClass("summary.lme",
         representation(call = "call",
                        logLik = "logLik",
                        re = "summary.ssclme",
                        residuals = "numeric"))

setClass("summary.GLMM", representation(family = "family"),
         contains = "summary.lme")

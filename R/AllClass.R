## Class definitions for the package

setOldClass("data.frame")
setOldClass("family")
setOldClass("logLik")
setOldClass("terms")

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

setClass("VarCorr",
         representation(scale="numeric",
                        reSumry="list",
                        useScale="logical"),
         prototype = list(scale = 1.0, useScale = TRUE))

## Representation of a linear mixed effects model
setClass("lmer",
         representation(
                        flist = "list", # list of grouping factors
                        perm = "list",  # list of permutations of levels (0-based)
                        Parent = "list",# list of Parent arrays for ZZpO
                        D = "list",     # list of diagonal factors (upper triangle)
                        bVar = "list",  # list of conditional variance factors (upper triangle)
                        L = "list",     # list of blocks of L
                        ZZpO = "list",  # list of diagonal blocks of Z'Z+Omega
                        Omega = "list", # list of relative precision matrices
                        REML = "logical", 
                        RXX = "matrix", # Augmented RXX component or its inverse
                        RZX = "matrix", # Augmented RZX component or its inverse
                        XtX = "matrix", # Original X'X matrix
                        ZtZ = "list",   # list of blocks of Z'Z
                        ZtX = "matrix", # Original Z'X matrix
                        cnames = "list",# column names of model matrices
                        devComp = "numeric", # Components of deviance
                        deviance = "numeric", # Current deviance (ML and REML)
                        nc = "integer", # number of columns in (augmented)
                                        # model matrices and number of observations
                        Gp = "integer", # Pointers to groups of rows in RZX
                        status = "logical",
                        call = "call",
                        terms = "terms",
                        assign = "integer",
                        fitted = "numeric",
                        residuals = "numeric",
                        frame = "data.frame"
                        ),
         validity = function(object) {
             .Call("lmer_validate", object, PACKAGE = "Matrix")
         })

setClass("summary.lmer",
         representation(useScale="logical",
                        showCorrelation="logical"),
         contains = "lmer")

setClass("lmList.confint", contains = "array")

setClass("lmer.ranef", representation(varFac = "list", stdErr =
                                      "numeric"),
         contains = "list")

setClass("lmer.ranef.confint", contains = "list")

setClass("lmer.coef", representation(varFac = "list", stdErr =
                                      "numeric"),
         contains = "list")

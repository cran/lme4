# Ensure that the stats4 and lattice packages are available
.onLoad <- function(lib, pkg)
{
    require("methods")
    require("stats4")
    require("lattice")
}

setClass("pdMat",      # parameterized positive-definite matrices
         representation(form="formula",    # a model-matrix formula
                        Names="character", # column (and row) names
                        param="numeric",   # parameter vector
                        Ncol="integer",    # number of columns
                        factor="matrix",   # factor of the pos-def matrix
                        logDet="numeric"   # logarithm of the absolute value
                        ## of the determinant of the factor (i.e. half
                        ## the logarithm of the determinant of the matrix)
                        ),
         prototype(form=formula(NULL),
                   Names=character(0),
                   param=numeric(0),
                   Ncol=as.integer(0),
                   factor=matrix(numeric(0),0,0),
                   logDet=numeric(0))
         )

#setClass("pdSymm", contains="pdMat")    # general symmetric pd matrices

#setClass("pdScalar", contains="pdMat") # special case of positive scalars
setClass("pdLogChol", contains="pdMat") # default parameterization
setClass("pdNatural", contains="pdMat") # log sd and logistic of correlation
#setClass("pdMatrixLog", contains="pdSymm") # matrix logarithm parameterization

setClass("pdDiag", contains="pdMat")    # diagonal pd matrices

setClass("pdIdent", contains="pdMat")   # positive multiple of the identity

setClass("pdCompSymm", contains="pdMat") # compound symmetric pd matrices

setClass("pdBlocked",                   # block-diagonal pd matrices
         representation("pdMat", components = "list"))

                       # positive-definite symmetric matrices as matrices
setClass("pdmatrix", contains="matrix")

                       # factors of positive-definite symmetric matrices
setClass("pdfactor", representation("matrix", logDet = "numeric"))

                       # correlation matrices and standard deviations
setClass("corrmatrix", representation("matrix", stdDev = "numeric"))

                       # a single level in the random effects structure
setClass("lmeLevel",
         representation(precision="pdMat", # the relative precision matrix
                        groups="factor",   # grouping factor for the level
                        columns="integer", # columns in model-matrix for level
                        parsInd="integer", # indices into the parameter vector
                        originalRows="list",
                        decomposedRows="list",
                        storedRows="list",
                        nrow="integer",
                        updateFactor="matrix", # used for EM update
                                               # and gradient calculation
                        hessianArray="array",
                        nlev="integer"))

                       # basic LME object representation
setClass("reStruct",
         representation(random="list", fixed="formula",
                        offset="numeric",
                        dirtyDecomposed="logical", useWeighted="logical",
                        dirtyStored="logical", dirtyBbetas="logical",
                        logLik="numeric",
                        analyticHessian="logical",
                        REML="logical",
                        reverseOrder="integer",
                        origOrder="integer",
                        original="matrix",
                        weighted="matrix", stored="matrix",
                        decomposed="matrix", bbetas="numeric",
                        dontCopy = "logical", assign.X = "ANY"),
         prototype=list(fixed = formula(NULL), dirtyBbetas = TRUE,
                        dirtyDecomposed=TRUE, REML=FALSE, dirtyStored=TRUE,
                        useWeighted=FALSE, logLik=as.numeric(NA),
                        dontCopy = FALSE, analyticHessian=FALSE))

setClass("lmeLevelList", contains="list")

setClass("summary.pdMat", representation(cor = "corrmatrix",
                                         structName = "character",
                                         noCorrelation = "logical",
                                         formula = "formula"),
         prototype=list(structName="", formula=formula(NULL)))

setClass("summary.reStruct",
         representation(fixed="formula",
                        coefficients="matrix",
                        scale="numeric",
                        denomDF="integer",
                        REML="logical",
                        ngrps="integer",
                        nobs="integer",
                        corFixed="corrmatrix",
                        reSumry="list",
                        useScale="logical",
                        showCorrelation="logical"
                        ))

setClass("summary.lme",
         representation(call = "call",
                        logLik = "logLik",
                        AIC = "numeric",
                        BIC = "numeric",
                        re ="summary.reStruct",
                        residuals = "numeric")
         )

setClass("summary.glmm", representation(method="character",
                                        family="character",
                                        link="character"),
         contains="summary.lme")

## Current versions of the methods package define

setClass("lme",
         representation(reStruct = "reStruct",
                        frame = "data.frame",
                        na.action = "ANY",
                        fitted = "numeric",
                        call = "call"),
         prototype(frame=data.frame(), fitted = numeric(0)))

## This is needed for the family slot of glmmStruct
setOldClass("family")

## Structure for fitting glmm classes
setClass("glmm",
         representation(family="family", # The glm family
                        origy="numeric",
                        n="numeric",
                        prior.weights="numeric",
                        init.weights="numeric",
                        init.y="numeric",
                        method="character"),
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










## deepayan experimenting with possible groupedData structures:



setClass("groupedData",
         representation(data = "data.frame",
                        formula = "formula",
                        outer = "formula",
                        inner = "formula",
                        labels = "list",
                        units = "list"))


## these are now obsolete: (nlme data are now in lme4 as S4 groupedData objects)

# ## Temporarily added so that groupedData objects are also data.frame objects
# setOldClass(c("nfnGroupedData", "nfGroupedData", "groupedData",
#               "data.frame"))

# setOldClass(c("nffGroupedData", "nfGroupedData", "groupedData",
#               "data.frame"))

# setOldClass(c("nmGroupedData", "groupedData", "data.frame"))


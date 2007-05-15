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
			flist = "list",    # list of grouping factors
			Zt = "dgCMatrix",  # sparse representation of Z'
			X = "matrix",	   # X
			y = "numeric",	   # y
			wts = "numeric",   # weights
                        ## do we need this for mer?
			wrkres = "numeric",# working residuals (copy of y for LMMs)
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
			status = "integer"
			),
         validity = function(object) .Call(mer_validate, object)
	)

## Representation of linear and generalized linear mixed effects model
setClass("lmer",
	 representation(frame = "data.frame",
                        call = "call",	   # call to model-fitting function
                        terms = "terms"),  # terms for fixed-effects
	 contains = "mer")

setClass("glmer",
	 representation(family = "family", # glm family
                        weights = "numeric"),
	 contains = "lmer")

## Representation of linear and generalized linear mixed effects model
setClass("lmer2",
	 representation(## original data
                        frame = "data.frame", # model frame or empty frame
                        call = "call",	    # matched call to model-fitting function
                        terms = "terms",    # terms for fixed-effects
			flist = "list",     # list of grouping factors
			ZXyt = "dgCMatrix", # sparse form of [Z;X;-y]'
			weights = "numeric",# can be of length 0 for constant wts
                        offset = "numeric", # can be of length 0 for 0 offset
			cnames = "list",    # column names of model matrices
			Gp = "integer",     # pointers to groups of rows in ZXyt
                        dims = "integer",   # dimensions and indicators
			## quantities that vary with Z, X, y, weights or offset
			A = "dsCMatrix",    # tcrossprod(ZXyt) (w. wts and offset)
			## slots that vary during the optimization
			ST = "list",        # list of TSST' rep of rel. var. mats
			L = "CHMfactor",    # sparse Cholesky factor of A*
			deviance = "numeric", # ML and REML deviance and components
			## Secondary slots only evaluated when requested.
			fixef = "numeric",
			ranef = "numeric"
			),
         validity = function(object) .Call(lmer2_validate, object)
         )

setClass("glmer2",
	 representation(family = "family",
                        X = "matrix",      # model matrix for fixed effects
                        eta = "numeric",   # linear predictor
                        mu = "numeric",    # inverse link of linear predictor
                        moff = "numeric",  # model offset, if any
                        pwts = "numeric",  # prior weights, if any
                        y = "numeric"),    # response
	 contains = "lmer2")

setClass("nlmer",
	 representation(## original data
                        frame = "data.frame", # model frame or empty frame
                        pnames = "character", # parameter names for nonlinear model
                        call = "call",	    # matched call to model-fitting function
                        terms = "terms",    # terms for fixed-effects
			flist = "list",     # list of grouping factors
                        Xt = "dgCMatrix",   # sparse form of X'
			Zt = "dgCMatrix",   # sparse form of Z'
                        y = "numeric",      # response
			weights = "numeric",# can be of length 0 for constant wts
			cnames = "list",    # column names of model matrices
			Gp = "integer",     # pointers to groups of columns in Z
                        dims = "integer",   # dimensions and indicators
			## quantities that vary with Z, X, y, weights or offset
			## slots that vary during the optimization
			ST = "list",        # list of TSST' rep of rel. var. mats
			L = "CHMfactor",    # sparse Cholesky factor of A*
			deviance = "numeric", # ML and REML deviance and components
			fixef = "numeric",
			ranef = "numeric"
			),
         validity = function(object) .Call(nlmer_validate, object)
         )


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

setClass("summary.lmer2", # the "lmer2" result ``enhanced'' :
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
	 contains = "lmer2")

setClass("summary.lmer", contains = c("summary.mer", "lmer"))

setClass("summary.glmer", contains = c("summary.mer", "glmer"))

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


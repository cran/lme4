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

setClass("mer",
	 representation(## original data
                        env = "environment",# evaluation env for nonlinear model
                        nlmodel = "call",# nonlinear model call
                        frame = "data.frame",# model frame (or empty frame)
                        call = "call",   # matched call
                        flist = "data.frame",  # list of grouping factors
                        X = "matrix",    # fixed effects model matrix
                        Zt = "dgCMatrix",# sparse form of Z'
                        pWt = "numeric",# prior weights,
                        offset = "numeric", # length 0 -> no offset
                        y = "numeric",   # response vector
###FIXME: Eliminate the cnames slot.  Put the names on the elements of the ST slot.
#                        cnames = "list", # row/column names of els of ST
                        Gp = "integer",  # pointers to row groups of Zt
                        dims = "integer",# dimensions and indicators
                        ## slots that vary during optimization
                        ST = "list", # list of TSST' rep of rel. var. mats
                        V = "matrix",    # gradient matrix
                        A = "dgCMatrix", # (ZTS)'
                        Cm = "dgCMatrix", # AH'G^{-1}W^{1/2} when s > 0
                        Cx = "numeric",  # x slot of Cm when s == 1 (full Cm not stored)
                        L = "CHMfactor", # Cholesky factor of weighted P(AA' + I)P'
                        deviance = "numeric", # ML and REML deviance and components
			fixef = "numeric",# fixed effects (length p)
			ranef = "numeric",# random effects (length q)
                        u = "numeric",   # orthogonal random effects (q)
                        eta = "numeric", # unbounded predictor
                        mu = "numeric",  # fitted values at current beta and b
                        muEta = "numeric",# d mu/d eta evaluated at current eta
                        var = "numeric", # conditional variances of Y
                        resid = "numeric",# raw residuals at current beta and b
                        sqrtXWt = "matrix",# sqrt of model matrix row weights
                        sqrtrWt = "numeric",# sqrt of weights used with residuals
                        RZX = "matrix", # dense sol. to L RZX = ST'ZtX = AX
                        RX = "matrix"),  # Cholesky factor of downdated X'X
         validity = function(object) .Call(mer_validate, object))

setClass("merMCMC",
         representation(
                        Gp = "integer",   # Gp slot from mer object
                        ST = "matrix",    # matrix of sampled ST pars
                        call = "call",    # matched call
                        deviance = "numeric",# vector of sampled deviances
                        dims = "integer", # dims from original mer object
                        fixef = "matrix", # matrix of sampled fixed effects pars
                        nc = "integer",   # number of columns per r.e. term
                        ranef = "matrix", # optional matrix of sampled r.e.
                        sigma = "matrix"  # sigma samples (may have 0 columns)
                        ),
         validity = function(object) .Call(merMCMC_validate, object))
                        
setClass("summary.mer",                 # Additional slots in a summary object
         representation(           
			methTitle = "character",
			logLik= "logLik",
			ngrps = "integer",
			sigma = "numeric", # scale, non-negative number
			coefs = "matrix",
			vcov = "dpoMatrix",
			REmat = "matrix",
			AICtab= "data.frame"),
         contains = "mer")

#setClass("ranef.mer", contains = "list")

#setClass("coef.mer", contains = "list")

setClass("sparseRasch", representation =
         list(dims = "integer",
              Zt = "dgCMatrix",
              y = "numeric",
              deviance = "numeric",
              offset = "numeric",
              L = "CHMfactor",
              fixef = "numeric",
              mu = "numeric",
              muEta = "numeric",
              pWt = "numeric",              
              resid = "numeric",
              sqrtrWt = "numeric",
              var = "numeric"),
         validity = function(object) TRUE)

setClass("merExt",
         representation(X0 = "matrix",    # original fixed effects model matrix
                        Zt0 = "dgCMatrix",# original sparse form of Z'
                        pars = "numeric", # additional parameters
                        y0 = "numeric"),  # original response vector
         contains = "mer"
         )

setClass("lmerStratVar",
         representation(sfac = "factor"),
         contains = "merExt")

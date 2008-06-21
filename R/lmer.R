# lmer, glmer and nlmer plus methods and utilities

if (0) {
### FIXME: Move this function to the stats package
rWishart <- function(n, df, invScal)
### Random sample from a Wishart distribution
    .Call(lme4_rWishart, n, df, invScal)
}

### Utilities for parsing the mixed model formula

findbars <- function(term)
### Return the pairs of expressions that separated by vertical bars
{
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    if (!is.call(term)) stop("term must be of class call")
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}

nobars <- function(term)
### Return the formula omitting the pairs of expressions that are
### separated by vertical bars
{
    if (!('|' %in% all.names(term))) return(term)
    if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
    if (length(term) == 2) {
	nb <- nobars(term[[2]])
	if (is.null(nb)) return(NULL)
	term[[2]] <- nb
	return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

subbars <- function(term)
### Substitute the '+' function for the '|' function
{
    if (is.name(term) || !is.language(term)) return(term)
    if (length(term) == 2) {
	term[[2]] <- subbars(term[[2]])
	return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name('|'))
	term[[1]] <- as.name('+')
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

subnms <- function(term, nlist)
### Substitute any names from nlist in term with 1
{
    if (!is.language(term)) return(term)
    if (is.name(term)) {
        if (any(unlist(lapply(nlist, get("=="), term)))) return(1)
        return(term)
    }
    stopifnot(length(term) >= 2)
    for (j in 2:length(term)) term[[j]] <- subnms(term[[j]], nlist)
    term
}

slashTerms <- function(x)
### Return the list of '/'-separated terms in an expression that
### contains slashes
{
    if (!("/" %in% all.names(x))) return(x)
    if (x[[1]] != as.name("/"))
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}

makeInteraction <- function(x)
### from a list of length 2 return recursive interaction terms
{
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}

expandSlash <- function(bb)
### expand any slashes in the grouping factors returned by findbars
{
    if (!is.list(bb)) return(expandSlash(list(bb)))
    ## I really do mean lapply(unlist(... - unlist returns a
    ## flattened list in this case
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
            return(lapply(unlist(makeInteraction(trms)),
                          function(trm) substitute(foo|bar,
                                                   list(foo = x[[2]],
                                                        bar = trm))))
        x
    }))
}

### Utilities used in lmer, glmer and nlmer

createCm <- function(A, s)
### Create the nonzero pattern for the sparse matrix Cm from A.
### ncol(A) is s * ncol(Cm).  The s groups of ncol(Cm) consecutive
### columns in A are overlaid to produce Cm.
{
    stopifnot(is(A, "dgCMatrix"))
    s <- as.integer(s)[1]
    if (s == 1L) return(A)
    if ((nc <- ncol(A)) %% s)
        stop(gettextf("ncol(A) = %d is not a multiple of s = %d",
                      nc, s))
    ncC <- as.integer(nc / s)
    TA <- as(A, "TsparseMatrix")
    as(new("dgTMatrix", Dim = c(nrow(A), ncC),
           i = TA@i, j = as.integer(TA@j %% ncC), x = TA@x),
       "CsparseMatrix")
}

### FIXME: somehow the environment of the mf formula does not have
### .globalEnv in its parent list.  example(Mmmec, package = "mlmRev")
### used to have a formula of ~ offset(log(expected)) + ... and the
### offset function was not found in eval(mf, parent.frame(2))
lmerFrames <- function(mc, formula, contrasts, vnms = character(0))
### Create the model frame, X, Y, wts, offset and terms

### mc - matched call of calling function
### formula - two-sided formula
### contrasts - contrasts argument
### vnms - names of variables to be included in the model frame
{
    mf <- mc
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]

    ## The model formula for evaluation of the model frame.  It looks
    ## like a linear model formula but includes any random effects
    ## terms and any names of parameters used in a nonlinear mixed model.
    frame.form <- subbars(formula)      # substitute `+' for `|'
    if (length(vnms) > 0)               # add the variables names for nlmer
        frame.form[[3]] <-
            substitute(foo + bar,
                       list(foo = parse(text = paste(vnms, collapse = ' + '))[[1]],
                            bar = frame.form[[3]]))

    ## The model formula for the fixed-effects terms only.
    fixed.form <- nobars(formula)       # remove any terms with `|'
    if (inherits(fixed.form, "name"))   # RHS is empty - use `y ~ 1'
        fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))

    ## attach the correct environment
    environment(fixed.form) <- environment(frame.form) <- environment(formula)

    ## evaluate a model frame
    mf$formula <- frame.form
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fe <- mf                            # save a copy of the call
    mf <- eval(mf, parent.frame(2))

    ## evaluate the terms for the fixed-effects only (used in anova)
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2)) # allow model.frame to update them

    ## response vector
    Y <- model.response(mf, "any")
    ## avoid problems with 1D arrays, but keep names
    if(length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if(!is.null(nm)) names(Y) <- nm
    }
    mt <- attr(fe, "terms")

    ## Extract X checking for a null model. This check shouldn't be
    ## needed because an empty formula is changed to ~ 1 but it can't hurt.
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts) else matrix(,NROW(Y),0)
    storage.mode(X) <- "double"      # when ncol(X) == 0, X is logical
    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL

    ## Extract the weights and offset.  For S4 classes we want the
    ## `not used' condition to be numeric(0) instead of NULL
    wts <- model.weights(mf); if (is.null(wts)) wts <- numeric(0)
    off <- model.offset(mf); if (is.null(off)) off <- numeric(0)

    ## check weights and offset
    if (any(wts < 0))
        stop(gettextf("negative weights not allowed"))
    if(length(off) && length(off) != NROW(Y))
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                      length(off), NROW(Y)))

    ## remove the terms attribute from mf
    attr(mf, "terms") <- mt
    list(Y = Y, X = X, wts = as.double(wts), off = as.double(off), mf = mf, fixef = fixef)
}

isNested <- function(f1, f2)
### Is f1 nested within f2?  That is, does every level of f1 occur
### in conjunction with exactly one level of f2?
{
    f1 <- as.factor(f1)
    f2 <- as.factor(f2)
    stopifnot(length(f1) == length(f2))
    sm <- as(new("ngTMatrix",
                 i = as.integer(f2) - 1L,
                 j = as.integer(f1) - 1L,
                 Dim = c(length(levels(f2)),
                 length(levels(f1)))),
             "CsparseMatrix")
    all(diff(sm@p) < 2)
}

lmerFactorList <- function(formula, mf, rmInt, drop)
### Create the list of grouping factors and the corresponding
### transposed model matrices.
### rmInt is a logical scalar indicating if the `(Intercept)` column
### should be removed before creating Zt
### drop is a logical scalar indicating if elements with numeric value
### 0 should be dropped from the sparse model matrices
{
    ## create factor list for the random effects
    bars <- expandSlash(findbars(formula[[3]]))
    if (!length(bars)) stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    fl <- lapply(bars,
                 function(x)
             {
                 ff <- eval(substitute(as.factor(fac)[,drop = TRUE],
                                       list(fac = x[[3]])), mf)
                 im <- as(ff, "sparseMatrix") # transpose of indicators
                 mm <- model.matrix(eval(substitute(~ expr, # model matrix
                                                    list(expr = x[[2]]))),
                                    mf)
                 if (rmInt) {
                     if (is.na(icol <- match("(Intercept)", colnames(mm)))) break
                     if (ncol(mm) < 2)
                         stop("lhs of a random-effects term cannot be an intercept only")
                     mm <- mm[ , -icol , drop = FALSE]
                 }
                 ans <- list(f = ff,
                             A = do.call(rBind,
                             lapply(seq_len(ncol(mm)), function(j) im)),
                             Zt = do.call(rBind,
                             lapply(seq_len(ncol(mm)),
                                    function(j) {im@x <- mm[,j]; im})),
                             ST = matrix(0, ncol(mm), ncol(mm),
                             dimnames = list(colnames(mm), colnames(mm))))
                 if (drop) {
                     ## This is only used for nlmer models.
                     ## Need to do something more complicated for A
                     ## here.  Essentially you need to create a copy
                     ## of im for each column of mm, im@x <- mm[,j],
                     ## create the appropriate number of copies,
                     ## prepend matrices of zeros, then rBind and drop0.
                     ans$A@x <- rep(0, length(ans$A@x))
                     ans$Zt <- drop0(ans$Zt)
                 }
                 ans
             })
    ## order terms by decreasing number of levels in the factor but don't
    ## change the order if this is already true
    nlev <- sapply(fl, function(el) length(levels(el$f)))
    if (any(diff(nlev)) > 0) fl <- fl[rev(order(nlev))]
    ## separate the terms from the factor list
    trms <- lapply(fl, "[", -1)
    names(trms) <- NULL
    fl <- lapply(fl, "[[", "f")
    attr(fl, "assign") <- seq_along(fl)
    ## check for repeated factors
    fnms <- names(fl)
    if (length(fnms) > length(ufn <- unique(fnms))) {
        ## check that the lengths of the number of levels coincide
        fl <- fl[match(ufn, fnms)]
        attr(fl, "assign") <- match(fnms, ufn)
    }
    names(fl) <- ufn
    list(trms = trms, fl = fl)
}

checkSTform <- function(ST, STnew)
### Check that the 'STnew' argument matches the form of ST.
{
    stopifnot(is.list(STnew), length(STnew) == length(ST),
              all.equal(names(ST), names(STnew)))
    lapply(seq_along(STnew), function (i)
           stopifnot(class(STnew[[i]]) == class(ST[[i]]),
                     all.equal(dim(STnew[[i]]), dim(ST[[i]]))))
    all(unlist(lapply(STnew, function(m) all(diag(m) > 0))))
}

lmerControl <- function(msVerbose = getOption("verbose"))
### Control parameters for lmer, glmer and nlmer
{
    list(
### FIXME: Should the user have control of maxIter and tolerance? If
### so, how should they be passed to the lmer_optimize C function?
         ## maxIter = as.integer(maxIter),
         ## tolerance = as.double(tolerance),
	 msVerbose = as.integer(msVerbose))# "integer" on purpose
}

VecFromNames <- function(nms, mode = "numeric")
### Generate a named vector of the given mode
{
    ans <- vector(mode = mode, length = length(nms))
    names(ans) <- nms
    ans[] <- NA
    ans
}

dimsNames <- c("nf", "n", "p", "q", "s", "np", "REML", "fTyp", "lTyp",
               "vTyp", "nest", "useSc", "nAGQ", "cvg")

devNames <- c("ML", "REML", "ldL2", "ldRX2", "sigmaML", "sigmaREML",
               "pwrss", "disc", "usqr", "wrss")

mkdims <- function(fr, FL, start, s = 1L)
### Create the standard versions of flist, Zt, Gp, ST, A, Cm,
### Cx, L and dd
{
    fl <- FL$fl
    asgn <- attr(fl, "assign")
    trms <- FL$trms
    ST <- lapply(trms, `[[`, "ST")
    Ztl <- lapply(trms, `[[`, "Zt")
    Zt <- do.call(rBind, Ztl)
    Zt@Dimnames <- vector("list", 2)
    Gp <- unname(c(0L, cumsum(sapply(Ztl, nrow))))
    .Call(mer_ST_initialize, ST, Gp, Zt)
    A <- do.call(rBind, lapply(trms, `[[`, "A"))
    rm(Ztl, FL)                         # because they could be large
    nc <- sapply(ST, ncol)         # of columns in els of ST
    Cm <- createCm(A, s)
    L <- .Call(mer_create_L, Cm)
    if (s < 2) Cm <- new("dgCMatrix")
    if (!is.null(start) && checkSTform(ST, start)) ST <- start

    ## record dimensions and algorithm settings
    dd <- VecFromNames(dimsNames, "integer")
    dd["nf"] <- length(ST)            # number of random-effects terms
    dd["n"] <- nrow(fr$mf)            # number of observations
    dd["p"] <- ncol(fr$X)             # number of fixed-effects coefficients
    dd["q"] <- nrow(Zt)               # number of random effects
    dd["s"] <- s
    nvc <- sapply(nc, function (qi) (qi * (qi + 1))/2) # no. of var. comp.
### FIXME: Check number of variance components versus number of
### levels in the factor for each term. Warn or stop as appropriate
    dd["np"] <- as.integer(sum(nvc))    # number of parameters in optimization
    dd["REML"] <- 0L                    # glmer and nlmer don't use REML
    dd["fTyp"] <- 2L                    # default family is "gaussian"
    dd["lTyp"] <- 5L                    # default link is "identity"
    dd["vTyp"] <- 1L                    # default variance function is "constant"
    ## check for nesting of factors
    dd["nest"] <- all(sapply(seq_along(fl)[-1],
                             function(i) isNested(fl[[i-1]], fl[[i]])))
    dd["useSc"] <- 1L                   # default is to use the scale parameter
    dd["nAGQ"] <- 1L                    # default is Laplace
    dd["cvg"]  <- 0L                    # no optimization yet attempted
    dev <- VecFromNames(devNames, "numeric")
    fl <- do.call(data.frame, c(fl, check.names = FALSE))
    attr(fl, "assign") <- asgn

    list(Gp = Gp, ST = ST, A = A, Cm = Cm, L = L, Zt = Zt,
         dd = dd, dev = dev, flist = fl)
}

famNms <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
            "poisson", "quasibinomial", "quasipoisson", "quasi")
linkNms <- c("logit", "probit", "cauchit", "cloglog", "identity",
             "log", "sqrt", "1/mu^2", "inverse")
varNms <- c("constant", "mu(1-mu)", "mu", "mu^2", "mu^3")

famType <- function(family)
{
    if (!(fTyp <- match(family$family, famNms, nomatch = 0)))
        stop(gettextf("unknown GLM family: %s",
                      sQuote(family$family), domain = "R-lme4"))
    if (!(lTyp <- match(family$link, linkNms, nomatch = 0)))
        stop(gettextf("unknown link: %s",
                      sQuote(family$link), domain = "R-lme4"))
    vNam <- switch(fTyp,
                   "mu(1-mu)",          # binomial
                   "constant",          # gaussian
                   "mu^2",              # Gamma
                   "mu^3",              # inverse.gaussian
                   "mu",                # poisson
                   "mu(1-mu)",          # quasibinomial
                   "mu",                # quasipoisson
                   family$varfun)       # quasi
    if (!(vTyp <- match(vNam, varNms, nomatch = 0)))
        stop(gettextf("unknown GLM family: %s",
                      sQuote(family$family), domain = "R-lme4"))
    c(fTyp = fTyp, lTyp = lTyp, vTyp = vTyp)
}

convergenceMessage <- function(cvg)
### Create the convergence message
{
    msg <- switch(as.character(cvg),
                  "3" = "X-convergence (3)",
                  "4" = "relative convergence (4)",
                  "5" = "both X-convergence and relative convergence (5)",
                  "6" = "absolute function convergence (6)",

                  "7" = "singular convergence (7)",
                  "8" = "false convergence (8)",
                  "9" = "function evaluation limit reached without convergence (9)",
                  "10" = "iteration limit reached without convergence (9)",
                  "14" = "storage has been allocated (?) (14)",

                  "15" = "LIV too small (15)",
                  "16" = "LV too small (16)",
                  "63" = "fn cannot be computed at initial par (63)",
                  "65" = "gr cannot be computed at initial par (65)")
    if (is.null(msg))
        msg <- paste("See PORT documentation.  Code (", cvg, ")", sep = "")
    msg
}

mer_finalize <- function(ans, verbose)
{
    .Call(mer_optimize, ans, verbose)
    if (ans@dims["cvg"] > 6) warning(convergenceMessage(ans@dims["cvg"]))
    .Call(mer_update_ranef, ans)
    .Call(mer_update_mu, ans)
    ans
}

## Modifications to lmer often involve modifying model matrices before
## creating and optimizing the mer object.  Everything past the model
## matrices is encapsulated in this function
lmer_finalize <- function(mc, fr, FL, start, REML, verbose)
{
    Y <- as.double(fr$Y)
    if (is.list(start) && all(sort(names(start)) == sort(names(FL))))
        start <- list(ST = start)
    if (is.numeric(start)) start <- list(STpars = start)
    dm <- mkdims(fr, FL, start[["ST"]])
    stopifnot(length(levels(dm$flist[[1]])) < length(Y))
### FIXME: A kinder, gentler error message may be in order.
### This checks that the number of levels in a grouping factor < n
### Only need to check the first factor because it is the one with
### the most levels.
    dm$dd["REML"] <- as.logical(REML)
    swts <- sqrt(unname(fr$wts))
    Cx <- numeric(0)
    if (length(swts))
        Cx <- (dm$A)@x
    p <- dm$dd["p"]
    n <- length(Y)

    ans <- new(Class = "mer",
               env = new.env(),
               nlmodel = (~I(x))[[2]],
               frame = fr$mf,
               call = mc,
               flist = dm$flist,
               X = fr$X,
               Zt = dm$Zt,
               pWt = unname(fr$wts),
               offset = unname(fr$off),
### FIXME: Should y retain its names? As it stands any row names in the
### frame are dropped.  Really?  Are they part of the frame slot (if not
### reduced to 0 rows)?
               y = unname(Y),
               Gp = unname(dm$Gp),
               dims = dm$dd,
               ST = dm$ST,
               A = dm$A,
               Cm = dm$Cm,
               Cx = Cx,
               L = dm$L,
               deviance = dm$dev,
               fixef = fr$fixef,
               ranef = numeric(dm$dd["q"]),
               u = numeric(dm$dd["q"]),
               eta = numeric(n),
               mu = numeric(n),
               resid = numeric(n),
               sqrtrWt = swts,
               sqrtXWt = as.matrix(swts),
               RZX = matrix(0, dm$dd["q"], p),
               RX = matrix(0, p, p))
    if (!is.null(stp <- start$STpars) && is.numeric(stp)) {
        STp <- .Call(mer_ST_getPars, ans)
        if (length(STp) == length(stp))
            .Call(mer_ST_setPars, ans, stp)
    }
    mer_finalize(ans, verbose)
}

### The main event
lmer <-
    function(formula, data, family = NULL, REML = TRUE,
             control = list(), start = NULL, verbose = FALSE,
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
### Linear Mixed-Effects in R
{
    mc <- match.call()
    if (!is.null(family)) {             # call glmer
        mc[[1]] <- as.name("glmer")
        return(eval(mc))
    }
    stopifnot(length(formula <- as.formula(formula)) == 3)

    fr <- lmerFrames(mc, formula, contrasts) # model frame, X, etc.
    FL <- lmerFactorList(formula, fr$mf, 0L, 0L) # flist, Zt
    if (!is.null(method <- list(...)$method)) {
        warning(paste("Argument", sQuote("methood"),
                      "is deprecated.  Use", sQuote("REML"),
                      "instead"))
        REML <- match.arg(method, c("REML", "ML")) == "REML"
    }
    lmer_finalize(mc, fr, FL, start, REML, verbose)
}

## for backward compatibility
lmer2 <-
    function(formula, data, family = NULL, REML = TRUE,
             control = list(), start = NULL, verbose = FALSE,
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
{
    .Deprecated("lmer")
    mc <- match.call()
    mc[[1]] <- as.name("lmer")
    eval.parent(mc)
}

glmer <-
function(formula, data, family = gaussian, start = NULL,
         verbose = FALSE, nAGQ = 1, subset, weights,
         na.action, offset, contrasts = NULL, model = TRUE,
         control = list(), ...)
### Fit a generalized linear mixed model
{
    mc <- match.call()
                                        # Evaluate and check the family
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    if(is.null(family$family)) stop("'family' not recognized")
    if(family$family == "gaussian" && family$link == "identity") {
        mc[[1]] <- as.name("lmer")      # use lmer not glmer
        mc$family <- NULL
        return(eval.parent(mc))
    }
    stopifnot(length(formula <- as.formula(formula)) == 3)

    fr <- lmerFrames(mc, formula, contrasts) # model frame, X, etc.
    offset <- wts <- NULL
    if (length(fr$wts)) wts <- fr$wts
    if (length(fr$off)) offset <- fr$off
    glmFit <- glm.fit(fr$X, fr$Y, weights = wts, # glm on fixed effects
                      offset = offset, family = family,
                      intercept = attr(attr(fr$mf, "terms"), "intercept") > 0)
    FL <- lmerFactorList(formula, fr$mf, 0L, 0L) # flist, Zt
    if (is.list(start) && all(sort(names(start)) == sort(names(FL))))
        start <- list(ST = start)
    if (is.numeric(start)) start <- list(STpars = start)
    dm <- mkdims(fr, FL, start[["ST"]])
    ft <- famType(glmFit$family)
    dm$dd[names(ft)] <- ft
    dm$dd["useSc"] <- as.integer(!(famNms[dm$dd["fTyp"] ] %in%
				   c("binomial", "poisson")))
    if (!is.null(method <- list(...)$method)) {
        msg <- paste("Argument", sQuote("method"),
                     "is deprecated.\nUse", sQuote("nAGQ"),
                     "to choose AGQ.  PQL is not available.")
        if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
            warning(msg)
        } else stop(msg)
    }
    if ((nAGQ <- as.integer(nAGQ)) < 1) nAGQ <- 1L
    dm$dd["nAGQ"] <- nAGQ
    y <- unname(as.double(glmFit$y))
    ##    dimnames(fr$X) <- NULL
    p <- dm$dd["p"]
    fixef <- fr$fixef
    fixef[] <- coef(glmFit)
    if (!is.null(ff <- start$fixef) && is.numeric(ff) &&
        length(ff) == length(fixef)) fixef <- ff

    ans <- new(Class = "mer",
               env = new.env(),
               nlmodel = (~I(x))[[2]],
               frame = if (model) fr$mf else fr$mf[0,],
               call = mc, flist = dm$flist,
               Zt = dm$Zt, X = fr$X, y = y,
               pWt = unname(glmFit$prior.weights),
               offset = unname(fr$off),
               Gp = unname(dm$Gp),
               dims = dm$dd, ST = dm$ST, A = dm$A,
               Cm = dm$Cm, Cx = (dm$A)@x, L = dm$L,
               deviance = dm$dev,
               fixef = fixef,
               ranef = numeric(dm$dd["q"]),
               u = numeric(dm$dd["q"]),
               eta = unname(glmFit$linear.predictors),
               mu = unname(glmFit$fitted.values),
               muEta = numeric(dm$dd["n"]),
               var = numeric(dm$dd["n"]),
               resid = unname(glmFit$residuals),
               sqrtXWt = as.matrix(numeric(dm$dd["n"])),
               sqrtrWt = numeric(dm$dd["n"]),
               RZX = matrix(0, dm$dd["q"], p),
               RX = matrix(0, p, p))
    if (!is.null(stp <- start$STpars) && is.numeric(stp)) {
        STp <- .Call(mer_ST_getPars, ans)
        if (length(STp) == length(stp))
            .Call(mer_ST_setPars, ans, stp)
    }
    cv <- do.call("lmerControl", control)
    if (missing(verbose)) verbose <- cv$msVerbose
    mer_finalize(ans, verbose)
    ans
}

nlmer <- function(formula, data, start = NULL, verbose = FALSE,
                  nAGQ = 1, subset, weights, na.action,
                  contrasts = NULL, model = TRUE, control = list(), ...)
### Fit a nonlinear mixed-effects model
{
    mc <- match.call()
    formula <- as.formula(formula)
    if (length(formula) < 3) stop("formula must be a 3-part formula")
    nlform <- as.formula(formula[[2]])
    if (length(nlform) < 3)
        stop("formula must be a 3-part formula")
    nlmod <- as.call(nlform[[3]])
    if (is.numeric(start)) start <- list(fixef = start)
    s <- length(pnames <- names(start$fixef))
    stopifnot(length(start$fixef) > 0, s > 0,
              inherits(data, "data.frame"), nrow(data) > 1)
### FIXME: Allow for the data argument to be missing.  What should the
### default be?
    if (any(pnames %in% names(data)))
        stop("parameter names must be distinct from names of the variables in data")
    anms <- all.vars(nlmod)
    if (!all(pnames %in% anms))
        stop("not all parameter names are used in the nonlinear model expression")

    if (!length(vnms <- setdiff(anms, pnames)))
        stop("there are no variables used in the nonlinear model expression")
    if ((nAGQ <- as.integer(nAGQ)) < 1) nAGQ <- 1L

    ## create a frame in which to evaluate the factor list
    fr <- lmerFrames(mc,
                     eval(substitute(foo ~ bar,
                                     list(foo = nlform[[2]],
                                          bar = subnms(formula[[3]],
                                          lapply(pnames, as.name))))),
                     contrasts, vnms)
    mf <- fr$mf
    env <- new.env()
    lapply(names(mf), function(nm) assign(nm, env = env, mf[[nm]]))
    n <- nrow(mf)
    lapply(pnames,
           function(nm) assign(nm, env = env, rep(start$fixef[[nm]],
                                   length.out = n)))

    n <- nrow(mf)
    mf <- mf[rep(seq_len(n), s), ]
    row.names(mf) <- seq_len(nrow(mf))
    ss <- rep.int(n, s)
    for (nm in pnames)
        mf[[nm]] <- rep.int(as.numeric(nm == pnames), ss)
                                        # factor list and model matrices
    FL <- lmerFactorList(substitute(foo ~ bar, list(foo = nlform[[2]],
                                                    bar = formula[[3]])),
                         mf, TRUE, TRUE)
    X <- as.matrix(mf[,pnames])
    rownames(X) <- NULL
    xnms <- colnames(fr$X)
    if (!is.na(icol <- match("(Intercept)",xnms))) xnms <- xnms[-icol]
### FIXME: The only times there would be additional columns in the
### fixed effects would be as interactions with parameter names and
### they must be constructed differently
#    if (length(xnms) > 0)
#        Xt <- cbind(Xt, fr$X[rep.int(seq_len(n), s), xnms, drop = FALSE])
    dm <- mkdims(fr, FL, start$STpars, s)
    p <- dm$dd["p"] <- length(start$fixef)
    n <- dm$dd["n"]
    if ((nAGQ <- as.integer(nAGQ)) < 1) nAGQ <- 1L
    dm$dd["nAGQ"] <- nAGQ

    ans <- new(Class = "mer",
               env = env,
               nlmodel = nlmod,
               frame = if (model) fr$mf else fr$mf[0,],
               call = mc,
               flist = dm$flist,
               X = X,
               Zt = dm$Zt,
               pWt = unname(sqrt(fr$wts)),
               offset = unname(fr$off),
               y = unname(as.double(fr$Y)),
               Gp = unname(dm$Gp),
               dims = dm$dd,
               ## slots that change during the iterations
               ST = dm$ST,
               V = matrix(0, n, s, dimnames = list(NULL, pnames)),
               A = dm$A,
               Cm = dm$Cm,
               L = dm$L,
               deviance = dm$dev,
               fixef = start$fixef,
               ranef = numeric(dm$dd["q"]),
               u = numeric(dm$dd["q"]),
               eta = numeric(n),
               mu = numeric(n),
               resid = numeric(n),
               sqrtXWt = matrix(0, n, s, dimnames = list(NULL, pnames)),
               sqrtrWt = unname(sqrt(fr$wts)),
               RZX = matrix(0, dm$dd["q"], p),
               RX = matrix(0, p, p)
               )
    .Call(mer_update_mu, ans)
### Add a check that the parameter names match the column names of gradient
    cv <- do.call("lmerControl", control)
    if (missing(verbose)) verbose <- cv$msVerbose
    mer_finalize(ans, verbose)
}

#### Extractors specific to mixed-effects models

setMethod("coef", signature(object = "mer"),
	  function(object, ...)
      {
          if (length(list(...)))
              warning(paste('arguments named "',
                            paste(names(list(...)), collapse = ", "),
                                  '" ignored', sep = ''))
          fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
          ref <- ranef(object)
          val <- lapply(ref, function(x) fef[rep(1, nrow(x)),,drop = FALSE])
          for (i in seq(a = val)) {
              refi <- ref[[i]]
              row.names(val[[i]]) <- row.names(refi)
              nmsi <- colnames(refi)
              if (!all(nmsi %in% names(fef)))
                  stop("unable to align random and fixed effects")
              for (nm in nmsi) val[[i]][[nm]] <- val[[i]][[nm]] + refi[,nm]
          }
          class(val) <- "coef.mer"
          val
#          new("coef.mer", val)
       })

setAs("mer", "dtCMatrix", function(from)
### Extract the L matrix
      as(from@L, "sparseMatrix"))

setMethod("fixef", signature(object = "mer"),
          function(object, ...)
### Extract the fixed effects
          object@fixef)

setMethod("ranef", signature(object = "mer"),
	  function(object, postVar = FALSE, drop = FALSE, ...)
### Extract the random effects
      {
          fl <- object@flist
          levs <- lapply(fl, levels)
          asgn <- attr(fl, "assign")
          Gp <- object@Gp
          ii <- lapply(diff(Gp), seq_len)
          rr <- object@ranef
          cn <- lapply(object@ST, colnames)
          ans <-
              lapply(split(lapply(seq_len(length(ii)),
                                  function(i)
                                  data.frame(matrix(rr[ii[[i]] + Gp[i]],
                                                    nc = length(cn[[i]]),
                                                    dimnames = list(
                                                    levs[[asgn[i]]], cn[[i]])),
                                             check.names = FALSE)),
                           asgn), function(lst) do.call(cbind, lst))
          names(ans) <- names(fl)
          class(ans) <- "ranef.mer"
          if (postVar) {
              if (length(fl) < length(Gp)) .NotYetImplemented
              pV <- .Call(mer_postVar, object)
              for (i in seq_along(ans))
                  attr(ans[[i]], "postVar") <- pV[[i]]
          }
          if (drop)
              ans <- lapply(ans, function(el)
                        {
                            if (ncol(el) > 1) return(el)
                            pv <- drop(attr(el, "postVar"))
                            el <- drop(as.matrix(el))
                            if (!is.null(pv))
                                attr(el, "postVar") <- pv
                            el
                        })
          ans
      })

print.ranef.mer <- function(x, ...) print(unclass(x), ...)
print.coef.mer <- function(x, ...) print(unclass(x), ...)

setMethod("sigma", signature(object = "mer"),
          function (object, ...) {
              dd <- object@dims
              if (!dd["useSc"]) return(1)
              object@deviance[if (dd["REML"]) "sigmaREML"
                              else "sigmaML"]
          })

setMethod("VarCorr", signature(x = "mer"),
	  function(x, ...)
### Create the VarCorr object of variances and covariances
      {
          sc <- sigma(x)
	  ans <- lapply(cc <- .Call(mer_ST_chol, x),
                        function(ch) {
                            val <- crossprod(sc * ch) # variance-covariance
                            stddev <- sqrt(diag(val))
                            correl <- t(val / stddev)/stddev
                            diag(correl) <- 1
                            attr(val, "stddev") <- stddev
                            attr(val, "correlation") <- correl
                            val
                        })
          fl <- x@flist
          names(ans) <- names(fl)[attr(fl, "assign")]
          attr(ans, "sc") <- if (x@dims["useSc"]) sc else NA
          ans
      })

#### Methods for standard extractors for fitted models

setMethod("anova", signature(object = "mer"),
	  function(object, ...)
      {
	  mCall <- match.call(expand.dots = TRUE)
	  dots <- list(...)
	  modp <- if (length(dots))
	      sapply(dots, is, "mer") | sapply(dots, is, "lm") else logical(0)
	  if (any(modp)) {		# multiple models - form table
	      opts <- dots[!modp]
	      mods <- c(list(object), dots[modp])
	      names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)],
				    as.character)
	      mods <- mods[order(sapply(lapply(mods, logLik, REML = FALSE),
					attr, "df"))]
	      calls <- lapply(mods, slot, "call")
	      data <- lapply(calls, "[[", "data")
	      if (any(data != data[[1]]))
		  stop("all models must be fit to the same data object")
	      header <- paste("Data:", data[[1]])
	      subset <- lapply(calls, "[[", "subset")
	      if (any(subset != subset[[1]]))
		  stop("all models must use the same subset")
	      if (!is.null(subset[[1]]))
		  header <-
		      c(header, paste("Subset", deparse(subset[[1]]), sep = ": "))
	      llks <- lapply(mods, logLik, REML = FALSE)
	      Df <- sapply(llks, attr, "df")
	      llk <- unlist(llks)
	      chisq <- 2 * pmax(0, c(NA, diff(llk)))
	      dfChisq <- c(NA, diff(Df))
	      val <- data.frame(Df = Df,
				AIC = sapply(llks, AIC),
				BIC = sapply(llks, BIC),
				logLik = llk,
				"Chisq" = chisq,
				"Chi Df" = dfChisq,
				"Pr(>Chisq)" = pchisq(chisq, dfChisq, lower = FALSE),
				row.names = names(mods), check.names = FALSE)
	      class(val) <- c("anova", class(val))
              attr(val, "heading") <-
                  c(header, "Models:",
                    paste(rep(names(mods), times = unlist(lapply(lapply(lapply(calls,
                                           "[[", "formula"), deparse), length))),
                         unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
                         sep = ": "))
	      return(val)
	  }
	  else { ## ------ single model ---------------------
              p <- object@dims["p"]
	      ss <- (object@RX[seq_len(p), p + 1L, drop = TRUE])^2
	      names(ss) <- names(object@fixef)
	      asgn <- attr(object@X, "assign")
	      terms <- terms(object)
	      nmeffects <- attr(terms, "term.labels")
	      if ("(Intercept)" %in% names(ss))
		  nmeffects <- c("(Intercept)", nmeffects)
	      ss <- unlist(lapply(split(ss, asgn), sum))
	      df <- unlist(lapply(split(asgn,  asgn), length))
	      #dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
	      ms <- ss/df
	      f <- ms/(sigma(object)^2)
	      #P <- pf(f, df, dfr, lower.tail = FALSE)
	      #table <- data.frame(df, ss, ms, dfr, f, P)
	      table <- data.frame(df, ss, ms, f)
	      dimnames(table) <-
		  list(nmeffects,
#			c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
		       c("Df", "Sum Sq", "Mean Sq", "F value"))
	      if ("(Intercept)" %in% nmeffects)
		  table <- table[-match("(Intercept)", nmeffects), ]
	      attr(table, "heading") <- "Analysis of Variance Table"
	      class(table) <- c("anova", "data.frame")
	      table
	  }
      })

if (FALSE) {
    setMethod("confint", signature(object = "mer"),
              function(object, parm, level = 0.95, ...)
              .NotYetImplemented()
              )
}

setMethod("deviance", signature(object="mer"),
	  function(object, REML = NULL, ...)
      {
          if (missing(REML) || is.null(REML) || is.na(REML[1]))
              REML <- object@dims["REML"]
          object@deviance[ifelse(REML, "REML", "ML")]
      })

setMethod("fitted", signature(object = "mer"),
          function(object, ...)
          napredict(attr(object@frame, "na.action"), object@mu))

setMethod("formula", signature(x = "mer"),
	  function(x, ...)
	  x@call$formula
	  )

setMethod("logLik", signature(object="mer"),
	  function(object, REML = NULL, ...)
### Extract the log-likelihood or restricted log-likelihood
      {
          dims <- object@dims
          if (is.null(REML) || is.na(REML[1]))
              REML <- dims["REML"]
          val <- -deviance(object, REML = REML)/2
          attr(val, "nall") <- attr(val, "nobs") <- dims["n"]
          attr(val, "df") <-
              dims["p"] + dims["np"] + as.logical(dims["useSc"])
          attr(val, "REML") <-  as.logical(REML)
          class(val) <- "logLik"
          val
      })

setMethod("residuals", signature(object = "mer"),
	  function(object, ...)
          napredict(attr(object@frame, "na.action"), object@resid))

setMethod("resid", signature(object = "mer"),
	  function(object, ...)
          napredict(attr(object@frame, "na.action"), object@resid))

setMethod("simulate", "mer",
          function(object, nsim = 1, seed = NULL, ...)
      {
	  if(!is.null(seed)) set.seed(seed)
	  if(!exists(".Random.seed", envir = .GlobalEnv))
	      runif(1)		     # initialize the RNG if necessary
          RNGstate <- .Random.seed
          dims <- object@dims
          etasim <- as.vector(object@X %*% fixef(object)) +  # fixed-effect contribution
              sigma(object) * (as(t(object@A) %*%    # random-effects contribution
                               matrix(rnorm(nsim * dims["q"]), nc = nsim),
                                  "matrix")
                               ## residual contribution
                               + matrix(rnorm(nsim * dims["n"]), nc = nsim))
          if (length(object@V) == 0 && length(object@muEta) == 0)
              return(etasim)
          stop("simulate method for GLMMs and NLMMs not yet implemented")
          })

setMethod("summary", signature(object = "mer"),
	  function(object, ...)
      {
          REML <- object@dims["REML"]
          fcoef <- fixef(object)
          vcov <- vcov(object)
          corF <- vcov@factors$correlation
          dims <- object@dims
          coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd) #, DF = DF)
          llik <- logLik(object, REML)
          dev <- object@deviance
          mType <- ifelse(non <- as.logical(length(object@V)), "NMM", "LMM")
          if (gen <- as.logical(length(object@muEta)))
              mType <- paste("G", mType, sep = '')
          mName <- switch(mType, LMM = "Linear", NMM = "Nonlinear",
                          GLMM = "Generalized linear",
                          GNMM = "Generalized nonlinear")
	  if(dims["nAGQ"] == 1)
              method <- "the Laplace approximation"
	  else
	      method <- "the adaptive Gaussian Hermite approximation"
          if (mType == "LMM")
              method <- ifelse(REML, "REML", "maximum likelihood")

          AICframe <- data.frame(AIC = AIC(llik), BIC = BIC(llik),
                                 logLik = c(llik),
                                 deviance = dev["ML"],
                                 REMLdev = dev["REML"],
                                 row.names = "")
          if (is.na(AICframe$REMLdev)) AICframe$REMLdev <- NULL
          varcor <- VarCorr(object)
          REmat <- formatVC(varcor)
          if (is.na(attr(varcor, "sc")))
              REmat <- REmat[-nrow(REmat), , drop = FALSE]

          if (nrow(coefs) > 0) {
              if (!dims["useSc"]) {
                  coefs <- coefs[, 1:2, drop = FALSE]
                  stat <- coefs[,1]/coefs[,2]
                  pval <- 2*pnorm(abs(stat), lower = FALSE)
                  coefs <- cbind(coefs, "z value" = stat, "Pr(>|z|)" = pval)
              } else {
                  stat <- coefs[,1]/coefs[,2]
                  ##pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
                  coefs <- cbind(coefs, "t value" = stat) #, "Pr(>|t|)" = pval)
              }
          } ## else : append columns to 0-row matrix ...
          new("summary.mer",
              object,
              methTitle = paste(mName, "mixed model fit by", method),
              logLik = llik,
              ngrps = sapply(object@flist, function(x) length(levels(x))),
              sigma = sigma(object),
              coefs = coefs,
              vcov = vcov,
              REmat = REmat,
              AICtab= AICframe
              )
      })## summary()

setMethod("model.frame", signature(formula = "mer"),
	  function(formula, ...) formula@frame)

setMethod("model.matrix", signature(object = "mer"),
	  function(object, ...) object@X)

setMethod("terms", signature(x = "mer"),
	  function(x, ...) attr(x@frame, "terms"))

setMethod("update", signature(object = "mer"),
	  function(object, formula., ..., evaluate = TRUE)
      {
	  call <- object@call
	  if (is.null(call))
	      stop("need an object with call slot")
	  extras <- match.call(expand.dots = FALSE)$...
	  if (!missing(formula.))
	      call$formula <- update.formula(formula(object), formula.)
	  if (length(extras) > 0) {
	      existing <- !is.na(match(names(extras), names(call)))
	      for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
	      if (any(!existing)) {
		  call <- c(as.list(call), extras[!existing])
		  call <- as.call(call)
	      }
	  }
	  if (evaluate)
	      eval(call, parent.frame())
	  else call
      })

setMethod("vcov", signature(object = "mer"),
	  function(object, ...)
### Extract the conditional variance-covariance matrix of the fixed effects
      {
          rr <- as(sigma(object)^2 *
                   chol2inv(object@RX, size = object@dims['p']), "dpoMatrix")
          nms <- colnames(object@X)
          dimnames(rr) <- list(nms, nms)
          rr@factors$correlation <- as(rr, "corMatrix")
          rr
      })

setMethod("with", signature(data = "mer"),
	  function(data, expr, ...) {
	      dat <- eval(data@call$data)
	      if (!is.null(na.act <- attr(data@frame, "na.action")))
		  dat <- dat[-na.act, ]
	      lst <- c(list(. = data), data@flist, data@frame, dat)
	      eval(substitute(expr), lst[unique(names(lst))])
	  })

### Show and print methods and utilities for them

formatVC <- function(varc, digits = max(3, getOption("digits") - 2))
### "format()" the 'VarCorr' matrix of the random effects -- for show()ing
{
    sc <- unname(attr(varc, "sc"))
    recorr <- lapply(varc, attr, "correlation")
    reStdDev <- c(lapply(varc, attr, "stddev"), list(Residual = sc))
    reLens <- unlist(c(lapply(reStdDev, length)))
    nr <- sum(reLens)
    reMat <- array('', c(nr, 4),
		   list(rep.int('', nr),
			c("Groups", "Name", "Variance", "Std.Dev.")))
    reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
    reMat[,2] <- c(unlist(lapply(varc, colnames)), "")
    reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
    reMat[,4] <- format(unlist(reStdDev), digits = digits)
    if (any(reLens > 1)) {
	maxlen <- max(reLens)
	corr <-
	    do.call("rBind",
		    lapply(recorr,
			   function(x, maxlen) {
			       x <- as(x, "matrix")
			       cc <- format(round(x, 3), nsmall = 3)
			       cc[!lower.tri(cc)] <- ""
			       nr <- dim(cc)[1]
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }, maxlen))
	colnames(corr) <- c("Corr", rep.int("", maxlen - 1))
	cbind(reMat, rBind(corr, rep.int("", ncol(corr))))
    } else reMat
}

## This is modeled a bit after  print.summary.lm :
printMer <- function(x, digits = max(3, getOption("digits") - 3),
                      correlation = TRUE, symbolic.cor = FALSE,
                      signif.stars = getOption("show.signif.stars"), ...)
{
    so <- summary(x)
    REML <- so@dims["REML"]
    llik <- so@logLik
    dev <- so@deviance
    dims <- x@dims

    cat(so@methTitle, "\n")
    if (!is.null(x@call$formula))
        cat("Formula:", deparse(x@call$formula),"\n")
    if (!is.null(x@call$data))
        cat("   Data:", deparse(x@call$data), "\n")
    if (!is.null(x@call$subset))
        cat(" Subset:",
            deparse(asOneSidedFormula(x@call$subset)[[2]]),"\n")
    print(so@AICtab, digits = digits)

    cat("Random effects:\n")
    print(so@REmat, quote = FALSE, digits = digits, ...)

    ngrps <- so@ngrps
    cat(sprintf("Number of obs: %d, groups: ", dims["n"]))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
    if (is.na(so@sigma))
	cat("\nEstimated scale (compare to 1):",
            sqrt(exp(so@deviance["lr2"])/so@dims["n"]), "\n")
    if (nrow(so@coefs) > 0) {
	cat("\nFixed effects:\n")
	printCoefmat(so@coefs, zap.ind = 3, #, tst.ind = 4
		     digits = digits, signif.stars = signif.stars)
	if(correlation) {
	    rn <- rownames(so@coefs)
	    corF <- so@vcov@factors$correlation
	    if (!is.null(corF)) {
		p <- ncol(corF)
		if (p > 1) {
		    cat("\nCorrelation of Fixed Effects:\n")
		    if (is.logical(symbolic.cor) && symbolic.cor) {
			print(symnum(as(corF, "matrix"), abbr.col = NULL))
		    }
		    else {
			corf <- matrix(format(round(corF@x, 3), nsmall = 3),
				       nc = p)
			dimnames(corf) <- list(abbreviate(rn, minlen=11),
					       abbreviate(rn, minlen=6))
			corf[!lower.tri(corf)] <- ""
			print(corf[-1, -p, drop=FALSE], quote = FALSE)
		    }
		}
	    }
	}
    }
    invisible(x)
}

setMethod("print", "mer", printMer)
setMethod("show", "mer", function(object) printMer(object))

printNlmer <- function(x, digits = max(3, getOption("digits") - 3),
                       correlation = TRUE, symbolic.cor = FALSE,
                       signif.stars = getOption("show.signif.stars"), ...)
### FIXME: Does nlmer need a separate show method?
{
    dims <- x@dims
    cat("Nonlinear mixed model fit by Laplace\n")
    if (!is.null(x@call$formula))
        cat("Formula:", deparse(x@call$formula),"\n")
    if (!is.null(x@call$data))
        cat("   Data:", deparse(x@call$data), "\n")
    if (!is.null(x@call$subset))
        cat(" Subset:",
            deparse(asOneSidedFormula(x@call$subset)[[2]]),"\n")

    cat("Random effects:\n")
    print(formatVC(VarCorr(x)), quote = FALSE,
          digits = max(3, getOption("digits") - 3))

    cat(sprintf("Number of obs: %d, groups: ", dims["n"]))
    ngrps <- sapply(x@flist, function(x) length(levels(x)))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
    cat("\nFixed effects:\n")
    print(x@fixef)
    invisible(x)
}

setMethod("refit", signature(object = "mer", newresp = "numeric"),
          function(object, newresp, ...)
      {
          newresp <- as.double(newresp[!is.na(newresp)])
          stopifnot(length(newresp) == object@dims["n"])
          object@y <- newresp
          mer_finalize(object, FALSE) # non-verbose fit
      })

setMethod("expand", signature(x = "mer"),
          function(x, sparse = TRUE, ...)
      {
          elexpand <- function(mat)
              list(T = new("dtrMatrix", uplo = "L", diag = "U",
                   x = as.vector(mat),
                   Dim = dim(mat), Dimnames = dimnames(mat)),
                   S = Diagonal(x = diag(mat)))
          if (!sparse) {
              ans <- x@ST
              fl <- x@flist
              if (all(attr(fl, "assign") == seq_along(ans)))
                  names(ans) <- names(fl)
              return(lapply(ans, elexpand))
          }
          .NotYetImplemented
      })

#### Methods for secondary, derived classes

setMethod("deviance", signature(object = "summary.mer"), function(object) object@deviance)
setMethod("logLik", signature(object = "summary.mer"), function(object) object@logLik)
setMethod("vcov", signature(object = "summary.mer"), function(object) object@vcov)
setMethod("summary", signature(object = "summary.mer"), function(object) object)

#### Methods to produce specific plots

plot.coef.mer <- function(x, y, ...)
{
    varying <- unique(do.call("c",
                              lapply(x, function(el)
                                     names(el)[sapply(el,
                                                      function(col)
                                                      any(col != col[1]))])))
    gf <- do.call("rBind", lapply(x, "[", j = varying))
    gf$.grp <- factor(rep(names(x), sapply(x, nrow)))
    switch(min(length(varying), 3),
           qqmath(eval(substitute(~ x | .grp,
                                  list(x = as.name(varying[1])))), gf, ...),
           xyplot(eval(substitute(y ~ x | .grp,
                                  list(y = as.name(varying[1]),
                                       x = as.name(varying[2])))), gf, ...),
           splom(~ gf | .grp, ...))
}

plot.ranef.mer <- function(x, y, ...)
{
    lapply(x, function(x) {
        cn <- lapply(colnames(x), as.name)
        switch(min(ncol(x), 3),
               qqmath(eval(substitute(~ x, list(x = cn[[1]]))), x, ...),
               xyplot(eval(substitute(y ~ x,
                                      list(y = cn[[1]],
                                           x = cn[[2]]))), x, ...),
               splom(~ x, ...))
    })
}

qqmath.ranef.mer <- function(x, data, ...)
{
    prepanel.ci <- function(x, y, se, subscripts, ...) {
        y <- as.numeric(y)
        se <- as.numeric(se[subscripts])
        hw <- 1.96 * se
        list(ylim = range(y - hw, y + hw, finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16, ...)  {
        panel.grid(h = -1,v = -1)
        panel.abline(h = 0)
        x <- as.numeric(x)
        y <- as.numeric(y)
        se <- as.numeric(se[subscripts])
        ly <- y - 1.96 * se
        uy <- y + 1.96 * se
        panel.segments(x, y - 1.96*se, x, y + 1.96 * se,
                       col = 'black')
        panel.xyplot(x, y, pch = pch, ...)
    }
    f <- function(x) {
        if (!is.null(pv <- attr(x, "postVar"))) {
            cols <- 1:(dim(pv)[1])
            se <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
            nr <- nrow(x)
            nc <- ncol(x)
            ord <- unlist(lapply(x, order)) +
                rep((0:(nc - 1)) * nr, each = nr)
            rr <- 1:nr
            ind <- gl(ncol(x), nrow(x), labels = names(x))
            xyplot(unlist(x)[ord] ~
                   rep(qnorm((rr - 0.5)/nr), ncol(x)) | ind[ord],
                   se = se[ord], prepanel = prepanel.ci, panel = panel.ci,
                   scales = list(y = list(relation = "free")),
                   xlab = "Standard normal quantiles",
                   ylab = NULL, aspect = 1, ...)
        } else {
            qqmath(~values|ind, stack(x),
                   scales = list(y = list(relation = "free")),
                   xlab = "Standard normal quantiles",
                   ylab = NULL, ...)
        }
    }
    lapply(x, f)
}


#### Creating and displaying a Markov Chain Monte Carlo sample from
#### the posterior distribution of the parameters

setMethod("mcmcsamp", signature(object = "mer"),
	  function(object, n = 1, verbose = FALSE, saveb = FALSE, ...)
### Generate a Markov chain Monte Carlo sample from the posterior distribution
### of the parameters in a linear mixed model
      {
          n <- max(1, as.integer(n)[1])
          dd <- object@dims
          ranef <- matrix(numeric(0), nrow = dd["q"], ncol = 0)
          if (saveb) ranef <- matrix(object@ranef, nrow = dd["q"], ncol = n)
          sigma <- matrix(unname(sigma(object)), nrow = 1,
                          ncol = (if (dd["useSc"]) n else 0))
          ff <- object@fixef
          fixef <- matrix(ff, dd["p"], n)
          rownames(fixef) <- names(ff)
          ans <- new("merMCMC",
                     Gp = object@Gp,
                     ST = matrix(.Call(mer_ST_getPars, object), dd["np"], n),
                     call = object@call,
                     dims = object@dims,
                     deviance = rep(unname(object@deviance["ML"]), n),
                     fixef = fixef,
                     nc = sapply(object@ST, nrow),
                     ranef = ranef,
                     sigma = sigma)
          .Call(mer_MCMCsamp, ans, object)
      })

setMethod("HPDinterval", signature(object = "merMCMC"),
          function(object, prob = 0.95, ...)
      {
          nms <- c("fixef", "ST")
          if (length(object@sigma)) nms <- c(nms, "sigma")
          if (length(object@ranef)) nms <- c(nms, "ranef")
          names(nms) <- nms
          lapply(lapply(nms, slot, object = object),
                 HPDinterval, prob = prob)
      })

setMethod("HPDinterval", signature(object = "matrix"),
          function(object, prob = 0.95, ...)
      {
          if (ncol(object) > nrow(object))
              object <- t(object)
          vals <- apply(object, 2, sort)
          if (!is.matrix(vals))
              stop("object must have nsamp > 1")
          nsamp <- nrow(vals)
          npar <- ncol(vals)
          gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
          init <- 1:(nsamp - gap)
          inds <- apply(vals[init + gap, , drop = FALSE] -
                        vals[init, , drop = FALSE], 2, which.min)
          ans <- cbind(vals[cbind(inds, 1:npar)],
                       vals[cbind(inds + gap, 1:npar)])
          dimnames(ans) <- list(colnames(object), c("lower", "upper"))
          attr(ans, "Probability") <- gap/nsamp
          ans
      })

### FIXME: Watch the names of the variance components here
setMethod("VarCorr", signature(x = "merMCMC"),
          function(x, type = c("raw", "varcov", "sdcorr", "logs"), ...)
      {
          if ("raw" == (type <- match.arg(type))) {
              ST <- t(x@ST)
              colnames(ST) <- paste("ST", 1:ncol(ST), sep = '')
              if (length(x@sigma)) return(cbind(ST, sigma = as.vector(x@sigma)))
              return(ST)
          }
          .Call(merMCMC_VarCorr, x, match(type, c("raw", "varcov", "sdcorr", "logs")))
      })

setMethod("as.matrix", signature(x = "merMCMC"),
          function(x, ...)
          cbind(t(x@fixef), VarCorr(x, ...)))

setMethod("as.data.frame", signature(x = "merMCMC"),
          function(x, row.names = NULL, optional = FALSE, ...)
          as.data.frame(as.matrix(x, ...), row.names = row.names, optional = optional, ...))

setAs("merMCMC", "data.frame", function(from) as.data.frame(from))

aslatticeframe <- function(x, ...)
{
    fr <- as.data.frame(x, ...)
    data.frame(dat = unlist(fr),
               par = gl(ncol(fr), nrow(fr), labels = colnames(fr)),
               iter = rep(1:nrow(fr), ncol(fr)))
}

## FIXME: More care should be taken to avoid duplicate argument names
## in the eventual call to lattice functions. Accumulate the arguments
## in a list and use do.call instead of direct calls.

setMethod("xyplot", signature(x = "merMCMC"),
          function(x, data, ...)
      {
          pfr <- aslatticeframe(x, ...)
          xyplot(dat ~ iter|par, pfr,
                 xlab = "Iteration number", ylab = NULL,
                 scales = list(x = list(axs = 'i'),
                 y = list(relation = "free", rot = 0)),
                 type = c("g", "l"),
                 layout = c(1, length(levels(pfr$par))),
                 strip = FALSE, strip.left = TRUE, ...)
      })

setMethod("densityplot", signature(x = "merMCMC"),
          function(x, data, ...)
          densityplot(~ dat | par, aslatticeframe(x, ...),
                      scales = list(relation = 'free'), ...)
          )

setMethod("qqmath", signature(x = "merMCMC"),
          function(x, data, ...)
          qqmath(~ dat | par, aslatticeframe(x, ...),
                 scales = list(y = list(relation = 'free')), ...)
          )


abbrvNms <- function(gnm, cnms)
### Abbreviate names of columns in grouping factors
### gnm - group name
### cnms - column names
{
    ans <- paste(abbreviate(gnm), abbreviate(cnms), sep = '.')
    if (length(cnms) > 1) {
	anms <- lapply(cnms, abbreviate, minlength = 3)
	nmmat <- outer(anms, anms, paste, sep = '.')
	ans <- c(ans, paste(abbreviate(gnm, minlength = 3),
			    nmmat[upper.tri(nmmat)], sep = '.'))
    }
    ans
}

mcmccompnames <- function(ans, object, saveb, trans, glmer, deviance)
### Mangle the names of the columns of the mcmcsamp result ans
### This operation is common to the methods for "lmer" and "glmer"
{
    gnms <- names(object@flist)
    cnms <- lapply(object@ST, colnames)
    ff <- fixef(object)
    colnms <- c(names(ff), if (glmer) character(0) else "sigma^2",
                unlist(lapply(seq_along(gnms),
                              function(i)
                              abbrvNms(gnms[i],cnms[[i]]))))
    if (trans) {
        ## parameter type: 0 => fixed effect, 1 => variance,
        ##		 2 => covariance
        ptyp <- c(integer(length(ff)), if (glmer) integer(0) else 1:1,
                  unlist(lapply(seq_along(gnms),
                                function(i)
                            {
                                k <- length(cnms[[i]])
                                rep(1:2, c(k, (k*(k-1))/2))
                            })))
        colnms[ptyp == 1] <-
            paste("log(", colnms[ptyp == 1], ")", sep = "")
        colnms[ptyp == 2] <-
            paste("atanh(", colnms[ptyp == 2], ")", sep = "")
    }
    if (deviance) colnms <- c(colnms, "deviance")
### FIXME: this will fail for a mer2 object
    if(saveb) {## maybe better colnames, "RE.1","RE.2", ... ?
        .NotYetImplemented()
        rZy <- object@rZy
        colnms <- c(colnms,
                    paste("b", sprintf(paste("%0",
                                             1+floor(log(length(rZy),10)),
                                             "d", sep = ''),
                                       seq_along(rZy)),
                          sep = '.'))
    }
    colnames(ans) <- colnms
    ans
}

devvals <- function(fm, pmat, sigma1 = FALSE)
{
    if (!is(fm, "mer"))
        stop('fm must be an "mer" fitted model')
### FIXME: add a check in here for glmer and nlmer
    np <- length(p0 <- .Call(mer_ST_getPars, fm))
    pmat <- as.matrix(pmat)
    if (ncol(pmat) != np + sigma1)
        stop(gettextf("pmat must have %d columns", np + sigma1))
    storage.mode(pmat) <- "double"
    if (is.null(pnms <- dimnames(pmat)[[2]]))
        pnms <- c(ifelse(sigma1, character(0), "sigma"),
                  paste("th", seq_len(np), sep = ""))
    dev <- fm@deviance
    ans <- matrix(0, nrow(pmat), ncol(pmat) + length(dev),
                  dimnames = list(NULL, c(pnms, names(dev))))
    for (i in seq_len(nrow(pmat))) {
        .Call(mer_ST_setPars, fm,
              ## This expression does not allow for correlated random
              ## effects.  It would be best to make the appropriate
              ## changes in the C code for ST_setPars.
              if (sigma1) pmat[i,-1]/pmat[i,1] else pmat[i,])
        ans[i, ] <- c(pmat[i, ], .Call(mer_update_RX, fm))
    }
    .Call(mer_ST_setPars, fm, p0)
    .Call(mer_update_RX, fm)
    as.data.frame(ans)
}

#### Odds and ends

## simulestimate <- function(x, FUN, nsim = 1, seed = NULL, control = list())
## {
##     FUN <- match.fun(FUN)
##     stopifnot((nsim <- as.integer(nsim[1])) > 0,
## 	      inherits(x, "lmer"))
##     if (!is.null(seed)) set.seed(seed)
##     ## simulate the linear predictors
##     lpred <- .Call(mer_simulate, x, nsim)
##     sc <- abs(x@devComp[8])
##     ## add fixed-effects contribution and per-observation noise term
##     lpred <- lpred + drop(x@X %*% fixef(x)) + rnorm(prod(dim(lpred)), sd = sc)

##     cv <- do.call(lmerControl, control)
##     Omega <- x@Omega
##     x@wrkres <- x@y <- lpred[,1]
##     .Call(mer_update_ZXy, x)
##     LMEoptimize(x) <- cv
##     template <- FUN(x)
##     if (!is.numeric(template))
##         stop("simulestimate currently only handles functions that return numeric vectors")
##     ans <- matrix(template, nr = nsim, nc = length(template), byrow = TRUE)
##     colnames(ans) <- names(template)
##     for (i in 1:nsim) {
##         x@wrkres <- x@y <- lpred[,i]
##         x@Omega <- Omega
##         .Call(mer_update_ZXy, x)
##         LMEoptimize(x) <- cv
##         foo <- try(FUN(x))
##         ans[i,] <- if (inherits(foo, "try-error")) NA else foo
##     }
##     ans
## }

hatTrace <- function(x)
{
    .NotYetImplemented()
    stopifnot(is(x, "mer"))
}

ST2Omega <- function(ST)
### Temporary function to convert the ST representation of the
### relative variance-covariance matrix returned by lmer into the
### Omega representation required by lmer
{
    if (nrow(ST) == 1) return(as(1/ST^2, "dpoMatrix"))
    dd <- diag(ST)
    T <- as(ST, "dtrMatrix")
    T@diag <- "U"
    crossprod(solve(T)/dd)
}


## setMethod("simulate", signature(object = "mer"),
## 	  function(object, nsim = 1, seed = NULL, ...)
##       {
## 	  if(!exists(".Random.seed", envir = .GlobalEnv))
## 	      runif(1)		     # initialize the RNG if necessary
## 	  if(is.null(seed))
## 	      RNGstate <- .Random.seed
## 	  else {
## 	      R.seed <- .Random.seed
## 	      set.seed(seed)
## 	      RNGstate <- structure(seed, kind = as.list(RNGkind()))
## 	      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
## 	  }

##           stopifnot((nsim <- as.integer(nsim[1])) > 0,
##                     inherits(object, "lmer"))
## 	  ## similate the linear predictors
## 	  lpred <- .Call(mer_simulate, object, nsim)
## 	  sc <- abs(object@devComp[8])

## 	  ## add fixed-effects contribution and per-observation noise term
## 	  lpred <- as.data.frame(lpred + drop(object@X %*% fixef(object)) +
## 				 rnorm(prod(dim(lpred)), sd = sc))
## 	  ## save the seed
## 	  attr(lpred, "seed") <- RNGstate
## 	  lpred
##       })

## We need to define an S4 print method, since using an S3 print
## method fails as soon as you call print() explicitly, e.g. when
## wanting to specify options.

## calculates degrees of freedom for fixed effects Wald tests
## This is a placeholder.  The answers are generally wrong.  It will
## be very tricky to decide what a 'right' answer should be with
## crossed random effects.

## setMethod("getFixDF", signature(object="mer"),
## 	  function(object, ...) {
## 	      devc <- object@devComp
## 	      rep(as.integer(devc[1]- devc[2]), devc[2])
## 	  })

## simss <- function(fm0, fma, nsim)
## {
##     ysim <- simulate(fm0, nsim)
##     cv <- list(gradient = FALSE, msMaxIter = 200:200,
## 	       msVerbose = 0:0)
##     sapply(ysim, function(yy) {
## 	.Call(mer_update_y, fm0, yy)
## 	LMEoptimize(fm0) <- cv
## 	.Call(mer_update_y, fma, yy)
## 	LMEoptimize(fma) <- cv
## 	exp(c(H0 = fm0@devComp[["logryy2"]],
## 	      Ha = fma@devComp[["logryy2"]]))
##     })
## }

## setMethod("denomDF", "mer",
##           function(x, ...)
##       {
##           mm <- x@X
##           aa <- attr(mm, "assign")
##           tt <- x@terms
##           if (!isNested(x))
##               return(list(coef = as.numeric(rep(NA, length(x@fixef))),
##                           terms = as.numeric(rep(NA,
##                           length(attr(tt, "order"))))))
##           hasintercept <- attr(tt, "intercept") > 0
##           ## check which variables vary within levels of grouping factors
##           vars <- eval(attr(tt, "variables"), x@frame)
##           fl <- x@flist
##           vv <- matrix(0:0, nrow = length(vars), ncol = length(fl),
##                         dimnames = list(NULL, names(fl)))
##           ## replace this loop by C code.
##           for (i in 1:nrow(ans))        # check if variables vary within factors
##               for (j in 1:ncol(ans))
##                   ans[i,j] <- all(tapply(vars[[i]], fl[[j]],
##                                          function(x) length(unique(x)) == 1))
##           ## which terms vary within levels of which grouping factors?
##           tv <- crossprod(attr(tt, "factors"), !ans)
##           ## maximum level at which the term is constant
##           ml <- apply(tv, 1, function(rr) max(0, which(as.logical(rr))))
##           ## unravel assignment applied to terms
##           ll <- attr(tt, "term.labels")
##           if (hasintercept)
##               ll <- c("(Intercept)", ll)
##           aaa <- factor(aa, labels = ll)
##           asgn <- split(order(aa), aaa)
##           nco <- lapply(asgn, length)   # number of coefficients per term
##           nlev <- lapply(fl, function(x) length(levels(x)))
##           if (hasintercept) asgn$"(Intercept)" <- NULL
##           list(ml = ml, nco = nco, nlev = nlev)
##       })

## Utilities for the fitted mer object
slotsz <- function(obj)
    rev(sort(sapply(slotNames(obj), function(s) object.size(slot(obj, s)))))

slotApply <- function(object, f, ..., simplify = FALSE) {
   .localFun <- function(what, ...) f(slot(object, what), ...)
   sapply(slotNames(object), .localFun, ..., simplify = simplify)
}


yfrm <- function(fm)
{
    stopifnot(is(fm, "mer"))
    snr <- slotApply(fm, function(x)
                 {
                     if (is(x, "matrix") ||
                         is(x, "data.frame") ||
                         is(x, "numeric")) return (NROW(x))
                     0
                 }, simplify = TRUE)
    snr <- snr[snr > 0 & !(names(snr) %in%
                           c("Gp", "dims", "deviance", "frame", "flist", "X"))]
    fr <- cbind(fm@frame, fm@flist[1:NROW(fm@frame), !(names(fm@flist) %in%
                                     names(fm@frame))])
    n <- NROW(fr)
    if (NROW(fm@X) == n)
        fr <- cbind(fr, X = fm@X, Xbeta = fm@X %*% fm@fixef,
                    Zb = crossprod(fm@Zt, fm@ranef)@x)
    do.call(cbind, c(list(fr), sapply(names(which(snr == NROW(fr))),
                                      slot, object = fm, simplify = FALSE)))
}

### Evaluate conditional components of a linear mixed model for a grid of ST
### parameter values.
###
### @param fm - a fitted linear mixed model
### @param parmat - a numeric matrix whose rows constitute suitable parameter
###     values for fm@ST
### @param type - which slot to extract
devmat <-
    function(fm, parmat, slotname = c("deviance", "fixef", "ranef", "u"), ...)
{
    stopifnot(is(fm, "mer"))
    dd <- fm@dims
    stopifnot(dd["fTyp"] == 2L, # gaussian family
              dd["lTyp"] == 5L, # identity link
              dd["vTyp"] == 1L, # variance function is "constant"
              length(fm@V) == 0L, # nonlinear parameter gradient is identity
              length(fm@muEta) == 0L) # eta -> mu map is identity
    oldpars <- .Call(mer_ST_getPars, fm)

    parmat <- as.matrix(parmat)
    storage.mode(parmat) <- "double"
    if (ncol(parmat) == dd["np"])
        parmat <- t(parmat)             # parameter vectors as columns
    stopifnot(nrow(parmat) == dd["np"])
    slotname <- match.arg(slotname)

    slotval <- function(x) {            # function to apply
        .Call(mer_ST_setPars, fm, x)
        .Call(mer_update_L, fm)
        .Call(mer_update_RX, fm)
        .Call(mer_update_ranef, fm)
        slot(fm, slotname)
    }
    ans <- apply(parmat, 2, slotval)
    slotval(oldpars)                    # restore the fitted model
    as.data.frame(t(rbind(parmat, ans)))
}

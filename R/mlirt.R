## Multilevel Iter Response Theory model

mlirt <-
    function(formula, data, family = binomial("probit"),
             control = list(), start = NULL,
             subset, weights, na.action = na.pass, offset,
             contrasts = NULL, model = FALSE, ...)
{
    formula <- as.formula(formula)
    if (length(formula) < 3) stop("formula must be a two-sided formula")
    cv <- do.call("lmerControl", control)

    ## Should difficulties be modeled as random effects?

    ranDiff <- ".item" %in% all.vars(formula)
    ## Establish model frame and fixed-effects model matrix and terms
    mc <- match.call()
    fr <- lmerFrames(mc, formula, data, contrasts)
    Y <- fr$Y

    ## check for a binary matrix response
    if (!is.matrix(Y) || !is.numeric(Y) ||
        any(!(unique(as.vector(Y)) %in% c(0, 1, NA)))) 
        stop("Response must be a binary, numeric matrix")
    nr <- nrow(Y)
    nc <- ncol(Y)
    ## expand model frame etc according to the items
    ind <- rep.int(1:nr, nc)
    mf <- fr$mf[ind, ]
    X <- fr$X[ind, ]
    if (ranDiff) mf$.item <- gl(nc, nr)
    else X <- cbind(X, -contr.sum(nc)[rep(1:nc,each = nr),])
    Y <- as.vector(Y)
    offset <- fr$offset[ind]

    mf$.subj <- gl(nr, 1, nr*nc)
    form <- substitute(Y ~ base + (1|.subj), list(base = formula[[3]])) 
    ## establish factor list and Ztl
    FL <- lmerFactorList(form, mf)
    cnames <- with(FL, c(lapply(Ztl, rownames), list(.fixed = colnames(X))))
    nc <- with(FL, sapply(Ztl, nrow))
    Ztl <- with(FL, .Call(Ztl_sparse, fl, Ztl))
    ## FIXME: change this when rbind has been fixed.
    Zt <- if (length(Ztl) == 1) Ztl[[1]] else do.call("rbind", Ztl)
    fl <- FL$fl

    ## check and evaluate the family argument
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if(is.function(family)) family <- family()
    if(is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (family$family != "binomial")
        stop("family must be a binomial family")
    fltype <- mkFltype(family)

    ## initial fit of a glm to the fixed-effects only.
    glmFit <- glm.fit(X, Y, weights = fl$weights[ind],
                      offset = offset, family = family,
                      intercept = attr(fr$mt, "intercept") > 0)
    Y <- as.double(glmFit$y)
    ## must do this after Y has possibly been reformulated
    mer <- .Call(mer_create, fl, Zt, X, Y, 0, nc, cnames)
    if (!is.null(start)) mer <- setOmega(mer, start)

    gVerb <- getOption("verbose")

    ## extract some of the components of glmFit
    ## weights could have changed
    weights <- glmFit$prior.weights
    eta <- glmFit$linear.predictors
    linkinv <- quote(family$linkinv(eta))
    mu.eta <- quote(family$mu.eta(eta))
    mu <- family$linkinv(eta)
    variance <- quote(family$variance(mu))
    dev.resids <- quote(family$dev.resids(Y, mu, weights))
    doLMEopt <- quote(LMEopt(x = mer, value = cv))
    mer@devComp[8] <- -1
    mer@status["glmm"] <- as.integer(2) # always use Laplace
    GSpt <- .Call(glmer_init, environment(), fltype)
    PQLpars <- c(coef(glmFit), .Call(mer_coef, mer, 2))
    fixInd <- seq(ncol(X))
    ## pars[fixInd] == beta, pars[-fixInd] == theta
    ## indicator of constrained parameters
    const <- c(rep(FALSE, length(fixInd)),
               unlist(lapply(mer@nc[seq(along = fl)],
                             function(k) 1:((k*(k+1))/2) <= k)
                      ))
    devLaplace <- function(pars) .Call(glmer_devLaplace, pars, GSpt)
    rel.tol <- abs(0.01/devLaplace(PQLpars))
    cat(paste("relative tolerance set to", rel.tol, "\n"))

    optimRes <- nlminb(PQLpars, devLaplace,
                       lower = ifelse(const, 5e-10, -Inf),
                       control = list(trace = cv$msVerbose,
                       iter.max = cv$msMaxIter,
                       rel.tol = rel.tol))
    .Call(glmer_finalize, GSpt)
    new("glmer",
        new("lmer", mer,
            frame =  data.frame(),
            terms = fr$mt, call = match.call()),
        weights = weights,
        family=family)
}

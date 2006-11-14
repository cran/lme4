## Preliminary version of a function to fit a model with carryover for
## the effect of one grouping factor (e.g. teacher) on another
## grouping factor (e.g. student)

carryOver <- function(formula, data, carry, REML = TRUE, control = list(),
                      start = NULL, varycarry = TRUE,
                      subset, weights, na.action, offset,
                      contrasts = NULL, model = TRUE, ...)
{
    formula <- as.formula(formula)
    if (length(formula) != 3) stop("formula must be a two-sided formula")
    cv <- do.call("lmerControl", control)

    ## Establish model frame and fixed-effects model matrix and terms
    mc <- match.call()
    fr <- lmerFrames(mc, formula, data, contrasts)
    Y <- fr$Y; X <- fr$X; weights <- fr$weights; offset <- fr$offset
    mf <- fr$mf; mt <- fr$mt
    
    ## establish factor list and Ztl
    FL <- lmerFactorList(formula, mf)
    fl <- FL$fl

    ## parse the carry-over formula
    carry <- as.formula(carry)
    if (length(carry) != 3) stop("carry must be a two-sided formula")
    if (!is.name(tvar <- carry[[2]]) ||
        !match(as.character(tvar), names(mf), nomatch = 0))
        stop("LHS of carry must be a name of a variable in formula")
    tvar <- eval(tvar, mf)
    if (!is.language(op <- carry[[3]]) || as.character(op[[1]]) != "/")
        stop("RHS of carry must be an expression of the form 'inner/outer'")
    if (!all(match(lapply(op[2:3], as.character), names(fl), nomatch = 0)))
        stop("Variables on RHS of carry must be names of grouping factors")
    outer <- eval(op[[3]], mf)
    if (is.factor(tvar)) tvar <- as.integer(tvar) - 1

    ## check the ordering
    ord <- order(outer, tvar)
    if (any(diff(ord) < 0)) {
        Y <- Y[ord]; X <- X[ord,]; mf <- mf[ord,]
        weights <- weights[ord]; offset <- offset[ord]
        FL <- lmerFactorList(formula, mf)
        fl <- FL$fl
        outer <- outer[ord]
        tvar <- tvar[ord]
    }
    nyr <- 1 + max(tapply(tvar, outer, function(x) diff(range(x))))
    disc0 <- 0.5^(0:(nyr - 1))
    disc <- disc0[-1]
    
    Ztsp <- .Call(Ztl_sparse, fl, FL$Ztl)
    innm <- as.character(op[[2]])
    Ztin <- Ztsp[[innm]]
    Ztsp[[innm]] <- .Call(Zt_carryOver, outer, Ztin, tvar, disc0)
    mer <- with(FL, .Call(mer_create, fl, do.call("rbind", Ztsp), X, Y, REML,
                          sapply(Ztl, nrow), # nc
                          c(lapply(Ztl, rownames), list(.fixed = colnames(X)))))
    if (!is.null(start)) mer <- setOmega(mer, start)
    if (cv$msMaxIter < 1) stop("msMaxIter must be positive")
    nc <- mer@nc
    np <- sum(unlist(lapply(nc, function(k) (k*(k+1))/2)))
    constr <- c(unlist(lapply(nc, function(k) 1:((k*(k+1))/2) <= k)),
                rep.int(FALSE, length(disc)))
    fn <- function(pars) {
        Ztsp[[innm]] <- .Call(Zt_carryOver, outer, Ztin, tvar,
                              c(1,pars[-(1:np)]))
        mer@Zt <- do.call("rbind", Ztsp)
        .Call(mer_update_ZXy, mer)
        deviance(.Call(mer_coefGets, mer, pars[1:np], 2))
    }
    start <- c(.Call(mer_coef, mer, 2), disc) #starting values
    fval <- fn(start)
    optimRes <- nlminb(start, fn, NULL,
                       lower = ifelse(constr, 5e-10, -Inf),
                       control = list(iter.max = cv$msMaxIter,
                       trace = as.integer(cv$msVerbose),
                       rel.tol = abs(0.001/fval)))
    estPar <- optimRes$par
    Ztsp[[innm]] <- .Call(Zt_carryOver, outer, Ztin, tvar,
                          c(1,estPar[-(1:np)]))
    mer@Zt <- do.call("rbind", Ztsp)
    .Call(mer_update_ZXy, mer)
    .Call(mer_coefGets, mer, estPar[1:np], 2)
    

    ## check for convergence on boundary
    if (any(bd <- (estPar[constr] < 1e-9))) {
        bpar <- rep.int(FALSE, length(estPar))
        bpar[constr] <- bd
        bgrp <- split(bpar,
                      rep(seq(along = nc),
                          unlist(lapply(nc,
                                        function(k) (k*(k+1))/2))))
        bdd <- unlist(lapply(bgrp, any))
        lens <- unlist(lapply(bgrp, length))
        if (all(lens[bdd] == 1)) { # variance components only
            warning("Estimated variance for ",
                    factorNames2char(names(mer@flist)[bdd]),
                    " is effectively zero\n")
        } else {
            warning("Estimated variance-covariance for ",
                    factorNames2char(names(mer@flist)[bdd]),
                    " is singular\n")
        }
    }
    if (optimRes$convergence != 0) {
        warning("nlminb returned message ", optimRes$message,"\n")
    }
    new("lmer", mer, frame = if (model) fr$mf else data.frame(),
        terms = mt, call = mc)
}

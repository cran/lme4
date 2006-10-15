
carryOver <- function(formula, data, carry, REML = TRUE, control = list(),
                      start = NULL, subset, weights, na.action, offset,
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
    FL <- lmerFactorList(formula, mf, X)
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

    ## check the ordering
    ord <- order(outer, tvar)
    if (any(diff(ord) < 0)) {
        warning("It is an advantage to have the data ordered by ",
                as.character(op[[3]]), " then ", as.character(carry[[2]]),
                " within ", as.character(op[[3]]), "\n")
        Y <- Y[ord]; X <- X[ord,]; mf <- mf[ord,]
        weights <- weights[ord]; offset <- offset[ord]
        FL <- lmerFactorList(formula, mf, X)
        fl <- FL$fl
        outer <- outer[ord]
    }

    Zt <- .Call(Zt_create_carryover, fl, FL$Ztl, 
                as.character(op[[2]]), as.character(op[[3]]))
    mer <- .Call(mer_create, fl, Zt, X, Y, REML, FL$nc, FL$cnames)
    if (!is.null(start)) mer <- setOmega(mer, start)
    .Call(mer_ECMEsteps, mer, cv$niterEM, cv$EMverbose)
    LMEoptimize(mer) <- cv
    new("lmer", mer, frame = if (model) fr$mf else data.frame(),
        terms = mt, call = mc)
}

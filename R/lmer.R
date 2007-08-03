# Methods for lmer and for the objects that it produces

## To Do: Determine why the names of the components of the values of
##   the ranef and coef extractor methods are not printed.

## Some utilities

## Return the pairs of expressions separated by vertical bars
findbars <- function(term)
{
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    if (!is.call(term)) stop("term must be of class call")
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}

## Return the formula omitting the pairs of expressions
## that are separated by vertical bars
nobars <- function(term)
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

## Substitute the '+' function for the '|' function
subbars <- function(term)
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

## Substitute any names from nlist in term with 1
subnms <- function(term, nlist)
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

## Return the list of '/'-separated terms in an expression that
## contains slashes
slashTerms <- function(x) {
    if (!("/" %in% all.names(x))) return(x)
    if (x[[1]] != as.name("/"))
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}

## from a list of length 2 return recursive interaction terms
makeInteraction <- function(x) {
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}


factorNames2char <- function(nms, collapse = ", ") {
    ## utility in messages / print etc:
    nms <- sQuote(nms)
    if(length(nms) == 1) paste("factor", nms)
    else paste("factors", paste(nms, collapse = collapse))
}

## expand any slashes in the grouping factors returned by findbars
expandSlash <- function(bb) {
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

abbrvNms <- function(gnm, cnms)
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

## Control parameters for lmer
lmerControl <-
  function(maxIter = 200, # used in ../src/lmer.c only
	   tolerance = sqrt(.Machine$double.eps),# ditto
	   msMaxIter = 200,
	   ## msTol = sqrt(.Machine$double.eps),
	   ## FIXME:  should be able to pass tolerances to nlminb()
	   msVerbose = getOption("verbose"),
	   niterEM = 15,
	   EMverbose = getOption("verbose"),
	   PQLmaxIt = 30,# FIXME: unused; PQL currently uses 'maxIter' instead
	   usePQL = FALSE,
	   gradient = TRUE,
	   Hessian = FALSE # unused _FIXME_
	   )
{
    list(maxIter = as.integer(maxIter),
	 tolerance = as.double(tolerance),
	 msMaxIter = as.integer(msMaxIter),
	 ## msTol = as.double(msTol),
	 msVerbose = as.integer(msVerbose),# "integer" on purpose
	 niterEM = as.integer(niterEM),
	 EMverbose = as.logical(EMverbose),
	 PQLmaxIt = as.integer(PQLmaxIt),
	 usePQL = as.logical(usePQL),
	 gradient = as.logical(gradient),
	 Hessian = as.logical(Hessian))
}

## check for predefined families
mkFltype <- function(family)
{
    fltype <- 0                         # not a predefined type
    if (family$family == "gaussian" && family$link == "identity") fltype <- -1
    if (family$family == "binomial" && family$link == "logit") fltype <- 1
    if (family$family == "binomial" && family$link == "probit") fltype <- 2
    if (family$family == "poisson" && family$link == "log") fltype <- 3
    as.integer(fltype)
}

rWishart <- function(n, df, invScal)
    .Call(lme4_rWishart, n, df, invScal)

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
          new("coef.lmer", val)
       })

setMethod("plot", signature(x = "coef.lmer"),
          function(x, y, ...)
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
      })

setMethod("plot", signature(x = "ranef.lmer"),
	  function(x, y, ...)
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
      })

setMethod("with", signature(data = "lmer"),
	  function(data, expr, ...) {
	      dat <- eval(data@call$data)
	      if (!is.null(na.act <- attr(data@frame, "na.action")))
		  dat <- dat[-na.act, ]
	      lst <- c(list(. = data), data@flist, data@frame, dat)
	      eval(substitute(expr), lst[unique(names(lst))])
	  })

setMethod("terms", signature(x = "lmer"),
	  function(x, ...) x@terms)

## Utility functions used in lmer and variations on lmer

## Check that the 'start' argument matches the form of the Omega
## slot in mer.  If so, install start as the Omega slot.
setOmega <- function(mer, start)
{
    Om <- mer@Omega
    if (!is.list(start) || length(start) != length(Om) ||
        !all.equal(names(Om), names(start)))
        stop(paste("start must be a list of length", length(Om),
                   "with names\n", paste(names(Om), collapse = ',')))
    for (i in seq(along = start)) {
        if (class(start[[i]]) != class(Om[[i]]) ||
            !all.equal(dim(start[[i]]), dim(Om[[i]])))
            stop(paste("start[[", i, "]] must be of class '", class(Om[[i]]),
                       "' and dimension ", paste(dim(Om[[i]]), collapse = ','),
                       sep = ''))
    }
    mer@Omega <- start
    mer
}

## Create model frame, terms, weights and offset from an lmer formula
##
## mc - matched call to parent function
## formula - two-sided formula
## contrasts - contrasts argument
## vnms - names of variables to be included in the model frame
lmerFrames <- function(mc, formula, contrasts, vnms = character(0))
{
    ## Must evaluate the model frame first and then fit the glm using
    ## that frame.  Otherwise missing values in the grouping factors
    ## cause inconsistent numbers of observations.
    mf <- mc
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    frame.form <- subbars(formula)      # substitute `+' for `|'
    if (length(vnms) > 0)
        frame.form[[3]] <-
            substitute(foo + bar,
                       list(foo = parse(text = paste(vnms, collapse = ' + '))[[1]],
                            bar = frame.form[[3]]))
    fixed.form <- nobars(formula)       # remove any terms with `|'
    if (inherits(fixed.form, "name"))  # RHS is empty - use a constant
        fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))
    environment(fixed.form) <- environment(frame.form) <- environment(formula)

    ## evaluate a model frame for fixed and random effects
    mf$formula <- frame.form
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fe <- mf
    mf <- eval(mf, parent.frame(2))

    ## get the terms for the fixed-effects only
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2))
    mt <- attr(fe, "terms")         # allow model.frame to update them
    ## response vector
    Y <- model.response(mf, "any")
    ## avoid problems with 1D arrays, but keep names
    if(length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if(!is.null(nm)) names(Y) <- nm
    }
    ## null model support
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts) else matrix(,NROW(Y),0)

    weights <- model.weights(mf)
    offset <- model.offset(mf)
    ## check weights and offset
    if( !is.null(weights) && any(weights < 0) )
        stop("negative weights not allowed")
    if(!is.null(offset) && length(offset) != NROW(Y))
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                      length(offset), NROW(Y)), domain = NA)
    list(Y = Y, X = X, weights = weights, offset = offset, mt = mt, mf = mf)
}

## Create the list of grouping factors and corresponding model
## matrices.  The model matrices are transposed in the Ztl list

## FIXME: Should the check on fltype be in this function or an auxillary check?
lmerFactorList <- function(formula, mf, fltype)
{
    ## create factor list for the random effects
    bars <- expandSlash(findbars(formula[[3]]))
    if (!length(bars)) stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    fl <- lapply(bars,
                 function(x)
                 eval(substitute(as.factor(fac)[,drop = TRUE],
                                 list(fac = x[[3]])), mf))
    ## order factor list by decreasing number of levels
    nlev <- sapply(fl, function(x) length(levels(x)))
    nobs <- length(fl[[1]])
    if(any(nlev == 0))
        stop("resulting factor(s) with 0 levels in random effects part:\n ",
             paste(sQuote(names(nlev[nlev == 0])), collapse=", "))
    ## Max # of levels allowed for a grouping factor.
    ## Binomial glmms can have nlev == nobs.
    maxlev <- nobs - !(fltype %in% 1:3)
    if (any(nlev > maxlev))
        stop("number of levels in grouping factor(s) ",
             paste(sQuote(names(nlev[nlev > maxlev])), collapse=", "),
             " is too large")
    if (any(diff(nlev) > 0)) {
        ord <- rev(order(nlev))
        bars <- bars[ord]
        fl <- fl[ord]
    }
    ## create list of transposed model matrices for random effects
    Ztl <- lapply(bars, function(x)
                  t(model.matrix(eval(substitute(~ expr,
                                                 list(expr = x[[2]]))),
                                 mf)))
    list(fl = fl, Ztl = Ztl)
}

lmer <- function(formula, data, family = gaussian,
                 method = c("REML", "ML", "PQL", "Laplace", "AGQ"),
                 control = list(), start = NULL,
                 subset, weights, na.action, offset, contrasts = NULL,
                 model = TRUE, ...)
{
    method <- match.arg(method)
    formula <- as.formula(formula)
    if (length(formula) < 3) stop("formula must be a two-sided formula")
    cv <- do.call("lmerControl", control)

    ## Establish model frame and fixed-effects model matrix and terms
    mc <- match.call()
    fr <- lmerFrames(mc, formula, contrasts)
    Y <- fr$Y; X <- fr$X; weights <- fr$weights; offset <- fr$offset
    mf <- fr$mf; mt <- fr$mt

    ## check and evaluate the family argument
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    if(is.null(family$family)) stop("'family' not recognized")
    fltype <- mkFltype(family)

    ## establish factor list and Ztl
    FL <- lmerFactorList(formula, mf, fltype)
    cnames <- with(FL, c(lapply(Ztl, rownames), list(.fixed = colnames(X))))
    nc <- with(FL, sapply(Ztl, nrow))
    Ztl <- with(FL, .Call(Ztl_sparse, fl, Ztl))
    ## FIXME: change this when rbind has been fixed.
    Zt <- if (length(Ztl) == 1) Ztl[[1]] else do.call("rBind", Ztl)
    fl <- FL$fl

    ## quick return for a linear mixed model
    if (fltype < 0) {
        mer <- .Call(mer_create, fl, Zt, X, as.double(Y),
                     method == "REML", nc, cnames)
        if (!is.null(start)) mer <- setOmega(mer, start)
        .Call(mer_ECMEsteps, mer, cv$niterEM, cv$EMverbose)
        LMEoptimize(mer) <- cv
        return(new("lmer", mer,
                   frame = if (model) fr$mf else data.frame(),
                   terms = mt,
                   call = mc))
    }

    ## The rest of the function applies to generalized linear mixed models
    if (method %in% c("ML", "REML")) method <- "Laplace"
    if (method == "AGQ")
        stop('method = "AGQ" not yet implemented for supernodal representation')
    if (method == "PQL") cv$usePQL <- TRUE # have to use PQL for method == "PQL"

    ## initial fit of a glm to the fixed-effects only.
    glmFit <- glm.fit(X, Y, weights = weights, offset = offset, family = family,
                      intercept = attr(mt, "intercept") > 0)
    Y <- as.double(glmFit$y)
    ## must do this after Y has possibly been reformulated
    mer <- .Call(mer_create, fl, Zt, X, Y, 0, nc, cnames)
    if (!is.null(start)) mer <- setOmega(mer, start)

    gVerb <- getOption("verbose")

    ## extract some of the components of glmFit
    ## weights could have changed
    weights <- glmFit$prior.weights
    if (is.null(offset)) offset <- numeric(length(Y))
    eta <- glmFit$linear.predictors
    linkinv <- quote(family$linkinv(eta))
    mu.eta <- quote(family$mu.eta(eta))
    mu <- family$linkinv(eta)
    variance <- quote(family$variance(mu))
    dev.resids <- quote(family$dev.resids(Y, mu, weights))
    LMEopt <- get("LMEoptimize<-")
    doLMEopt <- quote(LMEopt(x = mer, value = cv))
    if (family$family %in% c("binomial", "poisson")) # set the constant scale
        mer@devComp[8] <- -1
    mer@status["glmm"] <- as.integer(switch(method, PQL = 1, Laplace = 2, AGQ = 3))
    GSpt <- .Call(glmer_init, environment(), fltype)
    if (cv$usePQL) {
        .Call(glmer_PQL, GSpt)  # obtain PQL estimates
        PQLpars <- c(fixef(mer),
                     .Call(mer_coef, mer, 2))
    } else {
        PQLpars <- c(coef(glmFit),
                     .Call(mer_coef, mer, 2))
    }
    if (method == "PQL") {
        .Call(glmer_devLaplace, PQLpars, GSpt)
        .Call(glmer_finalize, GSpt)
        return(new("glmer",
                   new("lmer", mer,
                       frame = if (model) mf else data.frame(),
                       terms = mt, call = match.call()),
                   weights = weights,
                   family=family))
    }

    fixInd <- seq(ncol(X))
    ## pars[fixInd] == beta, pars[-fixInd] == theta
    ## indicator of constrained parameters
    const <- c(rep(FALSE, length(fixInd)),
               unlist(lapply(mer@nc[seq(along = fl)],
                             function(k) 1:((k*(k+1))/2) <= k)
                      ))
    devLaplace <- function(pars) .Call(glmer_devLaplace, pars, GSpt)
    rel.tol <- abs(0.01/devLaplace(PQLpars))
    if (cv$msVerbose) cat(paste("relative tolerance set to", rel.tol, "\n"))

    optimRes <- nlminb(PQLpars, devLaplace,
                       lower = ifelse(const, 5e-10, -Inf),
                       control = list(trace = cv$msVerbose,
                       iter.max = cv$msMaxIter,
                       rel.tol = rel.tol))
    .Call(glmer_finalize, GSpt)
    new("glmer",
        new("lmer", mer,
            frame = if (model) mf else data.frame(),
            terms = mt, call = match.call()),
        weights = weights,
        family=family)

}

## Extract the L matrix
setAs("mer", "dtCMatrix", function(from)
      .Call(mer_dtCMatrix, from))

## Extract the fixed effects
setMethod("fixef", signature(object = "mer"),
	  function(object, ...)
	  .Call(mer_fixef, object))

## Extract the random effects
setMethod("ranef", signature(object = "mer"),
	  function(object, postVar = FALSE, ...) {
	      ans <- new("ranef.lmer",
                         lapply(.Call(mer_ranef, object),
                                data.frame, check.names = FALSE))
              names(ans) <- names(object@flist)
              if (postVar) {
                  pV <- .Call(mer_postVar, object)
                  for (i in seq(along = ans))
                      attr(ans[[i]], "postVar") <- pV[[i]]
              }
              ans
	  })

## Optimization for mer objects
setReplaceMethod("LMEoptimize", signature(x="mer", value="list"),
		 function(x, value)
	     {
		 if (value$msMaxIter < 1) return(x)
		 nc <- x@nc
		 constr <- unlist(lapply(nc, function(k) 1:((k*(k+1))/2) <= k))
		 fn <- function(pars)
		     deviance(.Call(mer_coefGets, x, pars, 2))
                 start <- .Call(mer_coef, x, 2) #starting values
                 fval <- fn(start)
		 gr <- if (value$gradient)
		     function(pars) {
			 if (!isTRUE(all.equal(pars,
					       .Call(mer_coef, x,
						     2))))
			     .Call(mer_coefGets, x, pars, 2)
			 .Call(mer_gradient, x, 2)
		     }
		 else NULL
		 optimRes <- nlminb(.Call(mer_coef, x, 2),
				    fn, gr,
				    lower = ifelse(constr, 5e-10, -Inf),
				    control = list(iter.max = value$msMaxIter,
				    trace = as.integer(value$msVerbose),
                                    rel.tol = abs(0.001/fval)))
                 estPar <- optimRes$par
		 .Call(mer_coefGets, x, estPar, 2)

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
				 factorNames2char(names(x@flist)[bdd]),
				 " is effectively zero\n")
		     } else {
			 warning("Estimated variance-covariance for ",
				 factorNames2char(names(x@flist)[bdd]),
				 " is singular\n")
		     }
		 }
		 if (optimRes$convergence != 0) {
		     warning("nlminb returned message ", optimRes$message,"\n")
		 }
		 return(x)
	     })

setMethod("qqmath", signature(x = "ranef.lmer"),
          function(x, data, ...) {
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
          })

setMethod("deviance", signature(object = "mer"),
	  function(object, REML = NULL, ...) {
              if (is.null(REML)) REML <- object@status["REML"]
	      .Call(mer_factor, object)
	      object@deviance[[ifelse(REML, "REML", "ML")]]
	  })

## Mangle the names of the columns of the mcmcsamp result ans
## This operation is common to the methods for "lmer" and "glmer"
mcmccompnames <- function(ans, object, saveb, trans, glmer, deviance)
{
    gnms <- names(object@flist)
    cnms <- object@cnames
    ff <- fixef(object)
    colnms <- c(names(ff), if (glmer) character(0) else "sigma^2",
                unlist(lapply(seq(along = gnms),
                              function(i)
                              abbrvNms(gnms[i],cnms[[i]]))))
    if (trans) {
        ## parameter type: 0 => fixed effect, 1 => variance,
        ##		 2 => covariance
        ptyp <- c(integer(length(ff)), if (glmer) integer(0) else 1:1,
                  unlist(lapply(seq(along = gnms),
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
## FIXME: this will fail for a mer2 object
    if(saveb) {## maybe better colnames, "RE.1","RE.2", ... ?
        rZy <- object@rZy
        colnms <- c(colnms,
                    paste("b", sprintf(paste("%0",
                                             1+floor(log(length(rZy),10)),
                                             "d", sep = ''),
                                       seq(along = rZy)),
                          sep = '.'))
    }
    colnames(ans) <- colnms
    ans
}

setMethod("mcmcsamp", signature(object = "lmer"),
	  function(object, n = 1, verbose = FALSE, saveb = FALSE,
		   trans = TRUE, deviance = FALSE, ...)
      {
          ans <- t(.Call(mer_MCMCsamp, object, saveb, n, trans, verbose, deviance))
	  attr(ans, "mcpar") <- as.integer(c(1, n, 1))
	  class(ans) <- "mcmc"
	  mcmccompnames(ans, object, saveb, trans,
			glmer=FALSE, deviance=deviance)
      })

setMethod("mcmcsamp", signature(object = "glmer"),
	  function(object, n = 1, verbose = FALSE, saveb = FALSE,
		   trans = TRUE, deviance = FALSE, ...)
      {
          family <- object@family
          mer <- as(object, "mer")
          weights <- object@weights
          cv <- lmerControl()
          eta <- .Call(mer_fitted, mer)
          offset <- numeric(length(eta)) ## change this, save the offset in mer
          Y <- object@y
          linkinv <- quote(family$linkinv(eta))
          mu.eta <- quote(family$mu.eta(eta))
          mu <- family$linkinv(eta)
          variance <- quote(family$variance(mu))
          dev.resids <- quote(family$dev.resids(Y, mu, weights))
          LMEopt <- get("LMEoptimize<-")
          doLMEopt <- quote(LMEopt(x = mer, value = cv))
          fltype <- mkFltype(family)
          GSpt <- .Call(glmer_init, environment(), fltype)
          ans <- t(.Call(glmer_MCMCsamp, GSpt, saveb, n, trans, verbose, deviance))
          .Call(glmer_finalize, GSpt)
	  attr(ans, "mcpar") <- as.integer(c(1, n, 1))
	  class(ans) <- "mcmc"
	  mcmccompnames(ans, object, saveb, trans,
			glmer=TRUE, deviance=deviance)
      })

setMethod("simulate", signature(object = "mer"),
          function(object, nsim = 1, seed = NULL, ...)
      {
          if(!exists(".Random.seed", envir = .GlobalEnv))
              runif(1)                    # initialize the RNG if necessary
          if(is.null(seed))
              RNGstate <- .Random.seed
          else {
              R.seed <- .Random.seed
              set.seed(seed)
              RNGstate <- structure(seed, kind = as.list(RNGkind()))
              on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
          }
          
          stopifnot((nsim <- as.integer(nsim[1])) > 0, is(object, "lmer"))
          ## similate the linear predictors
          lpred <- .Call(lme4:::mer_simulate, object, nsim)
          sc <- abs(object@devComp[8])
          
          ## add fixed-effects contribution
          lpred <- lpred + drop(object@X %*% fixef(object))
          n <- prod(dim(lpred))
          
          if (!is(object, "glmer")) {  # and per-observation noise term
              response <- as.data.frame(lpred + rnorm(n, sd = sc))
              attr(response, "seed") <- RNGstate
              return(response)
          }
          fam <- attr(object, "family")
          if ("poisson"==fam$family) {
              response <- as.data.frame(matrix(byrow = FALSE, ncol=nsim,
                                          rpois(n, fam$linkinv(lpred))))
              attr(response, "seed") <- RNGstate
              return(response)
          }
          stop("simulate() not yet implemented for", fam$family,
               "glmms\n")
      })

### "FIXME": the following has too(?) much cut & paste from above

simulestimate <- function(x, FUN, nsim = 1, seed = NULL, control = list())
{
    FUN <- match.fun(FUN)
    stopifnot((nsim <- as.integer(nsim[1])) > 0,
              is(x, "lmer"))
    if (!is.null(seed)) set.seed(seed)
    ## simulate the linear predictors
    lpred <- .Call(mer_simulate, x, nsim)
    sc <- abs(x@devComp[8])
    ## add fixed-effects contribution and per-observation noise term
    lpred <- lpred + drop(x@X %*% fixef(x)) + rnorm(prod(dim(lpred)), sd = sc)

    cv <- do.call(lmerControl, control)
    Omega <- x@Omega
    x@wrkres <- x@y <- lpred[,1]
    .Call(mer_update_ZXy, x)
    LMEoptimize(x) <- cv
    template <- FUN(x)
    if (!is.numeric(template))
        stop("simulestimate currently only handles functions that return numeric vectors")
    ans <- matrix(template, nr = nsim, nc = length(template), byrow = TRUE)
    colnames(ans) <- names(template)
    for (i in 1:nsim) {
        x@wrkres <- x@y <- lpred[,i]
        x@Omega <- Omega
        .Call(mer_update_ZXy, x)
        LMEoptimize(x) <- cv
        foo <- try(FUN(x))
        ans[i,] <- if (inherits(foo, "try-error")) NA else foo
    }
    ans
}


formatVC <- function(varc, digits = max(3, getOption("digits") - 2))
{  ## "format()" the 'VarCorr'	matrix of the random effects -- for show()ing
    sc <- unname(attr(varc, "sc"))
    recorr <- lapply(varc, function(el) el@factors$correlation)
    reStdDev <- c(lapply(recorr, slot, "sd"), list(Residual = sc))
    reLens <- unlist(c(lapply(reStdDev, length)))
    nr <- sum(reLens)
    reMat <- array('', c(nr, 4),
		   list(rep.int('', nr),
			c("Groups", "Name", "Variance", "Std.Dev.")))
    reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
    reMat[,2] <- c(unlist(lapply(reStdDev, names)), "")
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

## We need to define an S4 print method, since using an S3 print
## method fails as soon as you call print() explicitly, e.g. when
## wanting to specify options.

## This is modeled a bit after  print.summary.lm :
printMer <- function(x, digits = max(3, getOption("digits") - 3),
                     correlation = TRUE, symbolic.cor = FALSE,
                     signif.stars = getOption("show.signif.stars"), ...)
{
    so <- summary(x)
    REML <- so@status["REML"]
    llik <- so@logLik
    dev <- so@deviance
    devc <- so@devComp
    glz <- so@isG

    cat(so@methTitle, "\n")
    if (!is.null(so@call$formula))
        cat("Formula:", deparse(so@call$formula),"\n")
    if (!is.null(so@call$data))
        cat("   Data:", deparse(so@call$data), "\n")
    if (!is.null(so@call$subset))
        cat(" Subset:",
            deparse(asOneSidedFormula(so@call$subset)[[2]]),"\n")
    if (glz)
        cat(" Family: ", so@family$family, "(",
            so@family$link, " link)\n", sep = "")
    print(so@AICtab, digits = digits)

    cat("Random effects:\n")
    print(so@REmat, quote = FALSE, digits = digits, ...)

    ngrps <- so@ngrps
    cat(sprintf("number of obs: %d, groups: ", devc[1]))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
    if (so@devComp[8] < 0)
	cat("\nEstimated scale (compare to ", abs(so@devComp[8]), ") ", so@sigma, "\n")
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
			corF <- matrix(format(round(corF@x, 3), nsmall = 3),
				       nc = p)
			dimnames(corF) <- list(abbreviate(rn, minlen=11),
					       abbreviate(rn, minlen=6))
			corF[!lower.tri(corF)] <- ""
			print(corF[-1, -p, drop=FALSE], quote = FALSE)
		    }
		}
	    }
	}
    }
    invisible(x)
}

setMethod("print", "mer", printMer)
setMethod("show", "mer", function(object) printMer(object))


setMethod("vcov", signature(object = "mer"),
	  function(object, REML = object@status["REML"], ...) {
	      rr <- as(object@devComp[8]^2 * tcrossprod(solve(object@RXX)), "dpoMatrix")
	      rr@factors$correlation <- as(rr, "corMatrix")
	      rr
	  })


## calculates degrees of freedom for fixed effects Wald tests
## This is a placeholder.  The answers are generally wrong.  It will
## be very tricky to decide what a 'right' answer should be with
## crossed random effects.

setMethod("getFixDF", signature(object="mer"),
	  function(object, ...) {
	      devc <- object@devComp
	      rep(as.integer(devc[1]- devc[2]), devc[2])
	  })

setMethod("logLik", signature(object="mer"),
	  function(object, REML = object@status["REML"], ...) {
	      val <- -deviance(object, REML = REML)/2
	      devc <- as.integer(object@devComp[1:2])
	      attr(val, "nall") <- attr(val, "nobs") <- devc[1]
	      attr(val, "df") <- abs(devc[2]) +
		  length(.Call(mer_coef, object, 0))
	      attr(val, "REML") <- REML
	      class(val) <- "logLik"
	      val
	  })

setMethod("VarCorr", signature(x = "mer"),
	  function(x, REML = x@status["REML"], ...)
      {
	  sc2 <- x@devComp[8]^2
	  cnames <- x@cnames
	  ans <- x@Omega
	  for (i in seq(a = ans)) {
	      el <- as(sc2 * solve(ans[[i]]), "dpoMatrix")
	      el@Dimnames <- list(cnames[[i]], cnames[[i]])
	      el@factors$correlation <- as(el, "corMatrix")
	      ans[[i]] <- el
	  }
	  attr(ans, "sc") <- sqrt(sc2)
	  ans
      })

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
				check.names = FALSE)
	      class(val) <- c("anova", class(val))
	      attr(val, "heading") <-
		  c(header, "Models:",
		    paste(names(mods),
			  unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
			  sep = ": "))
	      return(val)
	  }
	  else { ## ------ single model ---------------------
	      foo <- object
	      ss <- foo@rXy^2
	      ssr <- exp(foo@devComp["logryy2"])
	      names(ss) <- object@cnames[[".fixed"]]
	      asgn <- attr(foo@X, "assign")
	      terms <- foo@terms
	      nmeffects <- attr(terms, "term.labels")
	      if ("(Intercept)" %in% names(ss))
		  nmeffects <- c("(Intercept)", nmeffects)
	      ss <- unlist(lapply(split(ss, asgn), sum))
	      df <- unlist(lapply(split(asgn,  asgn), length))
	      #dfr <- getFixDF(object)
	      ms <- ss/df
	      #f <- ms/(ssr/dfr)
	      #P <- pf(f, df, dfr, lower.tail = FALSE)
	      #table <- data.frame(df, ss, ms, dfr, f, P)
	      table <- data.frame(df, ss, ms)
	      dimnames(table) <-
		  list(nmeffects,
#			c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
		       c("Df", "Sum Sq", "Mean Sq"))
	      if ("(Intercept)" %in% nmeffects)
		  table <- table[-match("(Intercept)", nmeffects), ]
	      attr(table, "heading") <- "Analysis of Variance Table"
	      class(table) <- c("anova", "data.frame")
	      table
	  }
      })

setMethod("confint", signature(object = "mer"),
	  function(object, parm, level = 0.95, ...)
	  .NotYetImplemented()
	  )

setMethod("fitted", signature(object = "mer"),
	  function(object, ...)
	  .Call(mer_fitted, object)
	  )

setMethod("formula", signature(x = "mer"),
	  function(x, ...)
	  x@call$formula
	  )

setMethod("residuals", signature(object = "glmer"),
	  function(object, ...) .NotYetImplemented())

setMethod("residuals", signature(object = "lmer"),
	  function(object, ...) object@y - fitted(object))

## FIXME: There should not be two identical methods like this but I'm not
##        sure how to pass the ... argument to a method for another generic
##        cleanly.
setMethod("resid", signature(object = "glmer"),
	  function(object, ...) .NotYetImplemented())

setMethod("resid", signature(object = "lmer"),
	  function(object, ...) object@y - fitted(object))

setMethod("summary", signature(object = "mer"),
	  function(object, ...) {

	      fcoef <- .Call(mer_fixef, object)
	      vcov <- vcov(object)
	      corF <- vcov@factors$correlation
	      ## DF <- getFixDF(object)
	      coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd) #, DF = DF)
	      REML <- object@status["REML"]
	      llik <- logLik(object, REML)
	      dev <- object@deviance
	      devc <- object@devComp

	      glz <- is(object, "glmer")
	      methTitle <-
		  if (glz)
		      paste("Generalized linear mixed model fit using",
			    switch(object@status["glmm"],
                                   "PQL", "Laplace", "AGQ"))
		  else paste("Linear mixed-effects model fit by",
			     if(REML) "REML" else "maximum likelihood")

	      AICframe <- {
		  if (glz)
		      data.frame(AIC = AIC(llik), BIC = BIC(llik),
				 logLik = c(llik),
				 deviance = -2*llik,
				 row.names = "")
		  else
		      data.frame(AIC = AIC(llik), BIC = BIC(llik),
				 logLik = c(llik),
				 MLdeviance = dev["ML"],
				 REMLdeviance = dev["REML"],
				 row.names = "")
	      }
	      REmat <- formatVC(VarCorr(object))
	      if (object@devComp[8] < 0) REmat <- REmat[-nrow(REmat), , drop = FALSE]

	      if (nrow(coefs) > 0) {
		  if (object@devComp[8] >= 0) {
		      stat <- coefs[,1]/coefs[,2]
		      ##pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
		      coefs <- cbind(coefs, "t value" = stat) #, "Pr(>|t|)" = pval)
		  } else {
		      coefs <- coefs[, 1:2, drop = FALSE]
		      stat <- coefs[,1]/coefs[,2]
		      pval <- 2*pnorm(abs(stat), lower = FALSE)
		      coefs <- cbind(coefs, "z value" = stat, "Pr(>|z|)" = pval)
		  }
	      } ## else : append columns to 0-row matrix ...

	      new(if(is(object, "glmer")) "summary.glmer" else
                  {if(is(object, "lmer")) "summary.lmer" else "summary.mer"},
		  object,
		  isG = glz,
		  methTitle = methTitle,
		  logLik = llik,
		  ngrps = sapply(object@flist, function(x) length(levels(x))),
		  sigma = .Call(mer_sigma, object, REML),
		  coefs = coefs,
		  vcov = vcov,
		  REmat = REmat,
		  AICtab= AICframe
		  )
	  })## summary()

setMethod("print", "mer", printMer)
setMethod("show", "mer", function(object) printMer(object))


## Methods for "summary.*" objects:
setMethod("vcov", signature(object = "summary.mer"),
	  function(object) object@vcov)
setMethod("logLik", signature(object = "summary.mer"),
	  function(object) object@logLik)
setMethod("deviance", signature(object = "summary.mer"),
 	  function(object) object@deviance)
setMethod("summary", signature(object = "summary.mer"), function(object) object)


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

simss <- function(fm0, fma, nsim)
{
    ysim <- simulate(fm0, nsim)
    cv <- list(gradient = FALSE, msMaxIter = 200:200,
	       msVerbose = 0:0)
    sapply(ysim, function(yy) {
	.Call(mer_update_y, fm0, yy)
	LMEoptimize(fm0) <- cv
	.Call(mer_update_y, fma, yy)
	LMEoptimize(fma) <- cv
	exp(c(H0 = fm0@devComp[["logryy2"]],
	      Ha = fma@devComp[["logryy2"]]))
    })
}

## Some leftover code from the old AGQ method in lmer.
if (FALSE) {
### FIXME: For nf == 1 change this to an AGQ evaluation.  Needs
### AGQ for nc > 1 first.
    fxd <- PQLpars[fixInd]
    loglik <- logLik(mer)

    if (method %in% c("Laplace", "AGQ")) {
	nAGQ <- 1
	if (method == "AGQ") {	  # determine nAGQ at PQL estimates
	    dev11 <- devAGQ(PQLpars, 11)
	    ## FIXME: Should this be an absolute or a relative tolerance?
	    devTol <- sqrt(.Machine$double.eps) * abs(dev11)
	    for (nAGQ in c(9, 7, 5, 3, 1))
		if (abs(dev11 - devAGQ(PQLpars, nAGQ - 2)) > devTol) break
	    nAGQ <- nAGQ + 2
	    if (gVerb)
		cat(paste("Using", nAGQ, "quadrature points per column\n"))
	}
	obj <- function(pars)
	    .Call(glmer_devAGQ, pars, GSpt, nAGQ)
	optimRes <-
	    nlminb(PQLpars, obj,
		   lower = ifelse(const, 5e-10, -Inf),
		   control = list(trace = getOption("verbose"),
		   iter.max = cv$msMaxIter))
	optpars <- optimRes$par
	if (optimRes$convergence != 0)
	    warning("nlminb failed to converge")
	deviance <- optimRes$objective
	if (gVerb)
	    cat(paste("convergence message", optimRes$message, "\n"))
	fxd[] <- optpars[fixInd]  ## preserve the names
	.Call(lmer_coefGets, mer, optpars[-fixInd], 2)
    }

    .Call(glmer_finalize, GSpt)
    loglik[] <- -deviance/2
}## end{leftover}

setMethod("isNested", "mer",
          function(x, ...) !(x@L@type[1]),
          valueClass = "logical")

setMethod("denomDF", "mer",
          function(x, ...)
      {
          mm <- x@X
          aa <- attr(mm, "assign")
          tt <- x@terms
          if (!isNested(x))
              return(list(coef = as.numeric(rep(NA, length(x@fixef))),
                          terms = as.numeric(rep(NA,
                          length(attr(tt, "order"))))))
          hasintercept <- attr(tt, "intercept") > 0
          ## check which variables vary within levels of grouping factors
          vars <- eval(attr(tt, "variables"), x@frame)
          fl <- x@flist
          vv <- matrix(0:0, nrow = length(vars), ncol = length(fl),
                        dimnames = list(NULL, names(fl)))
          ## replace this loop by C code.
          for (i in 1:nrow(ans))        # check if variables vary within factors
              for (j in 1:ncol(ans))
                  ans[i,j] <- all(tapply(vars[[i]], fl[[j]],
                                         function(x) length(unique(x)) == 1))
          ## which terms vary within levels of which grouping factors?
          tv <- crossprod(attr(tt, "factors"), !ans)
          ## maximum level at which the term is constant
          ml <- apply(tv, 1, function(rr) max(0, which(as.logical(rr))))
          ## unravel assignment applied to terms
          ll <- attr(tt, "term.labels")
          if (hasintercept)
              ll <- c("(Intercept)", ll)
          aaa <- factor(aa, labels = ll)
          asgn <- split(order(aa), aaa)
          nco <- lapply(asgn, length)   # number of coefficients per term
          nlev <- lapply(fl, function(x) length(levels(x)))
          if (hasintercept) asgn$"(Intercept)" <- NULL
          list(ml = ml, nco = nco, nlev = nlev)
      })

hatTrace <- function(x)
{
    stopifnot(is(x, "mer"))
    .Call(mer_hat_trace2, x)
}

## Check that the 'start' argument matches the form of the ST
## slot in mer2.  If so, install start as the ST slot.
setST <- function(mer, start)
{
    ST <- mer@ST
    if (!is.list(start) || length(start) != length(ST) ||
        !all.equal(names(ST), names(start)))
        stop(paste("start must be a list of length", length(ST),
                   "with names\n", paste(names(ST), collapse = ',')))
    for (i in seq(along = start)) {
        if (class(start[[i]]) != class(ST[[i]]) ||
            !all.equal(dim(start[[i]]), dim(ST[[i]])))
            stop(paste("start[[", i, "]] must be of class '", class(ST[[i]]),
                       "' and dimension ", paste(dim(ST[[i]]), collapse = ','),
                       sep = ''))
    }
    mer@ST <- start
    mer
}

lmer2 <- function(formula, data, family = gaussian,
                  method = c("REML", "ML", "Laplace", "AGQ"),
                  control = list(), start = NULL,
                  subset, weights, na.action, offset, contrasts = NULL,
                  model = TRUE, ...)
{
    method <- match.arg(method)
    formula <- as.formula(formula)
    if (length(formula) < 3) stop("formula must be a two-sided formula")
    cv <- do.call("lmerControl", control)

    ## Establish model frame and fixed-effects model matrix and terms
    mc <- match.call()
    fr <- lmerFrames(mc, formula, contrasts)
    storage.mode(fr$X) <- "double" # when ncol(X) == 0, X is logical

    ## Fit a glm to the fixed-effects only.
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    if(is.null(family$family)) stop("'family' not recognized")
    glmFit <- glm.fit(fr$X, fr$Y, weights = fr$weights,
                      offset = fr$offset, family = family,
                      intercept = attr(fr$mt, "intercept") > 0)
    storage.mode(glmFit$y) <- "double"

    ## establish factor list and Ztl
    FL <- lmerFactorList(formula, fr$mf, !is.null(family))
    Ztl <- with(FL, .Call(Ztl_sparse, fl, Ztl))

    mer <- .Call(lmer2_create, fr, FL, Ztl, glmFit,
                 method, mc, model)
    rm(fr, FL, Ztl, glmFit)
    if (!is.null(start)) mer <- setST(mer, start)
    if (!is(mer, "glmer2")) {      # linear mixed model
        .Call(lmer2_optimize, mer, cv$msVerbose)
        .Call(lmer2_update_effects, mer)
    }
    mer
}

## FIXME: The fixef and vcov methods should incorporate the permutation.
## Extract the fixed effects
setMethod("fixef", signature(object = "lmer2"),
          function(object, ...) {
              ans <- object@fixef
              names(ans) <- object@cnames$.fixed
              ans
          })

## Extract the conditional variance-covariance matrix of the fixed effects
setMethod("vcov", signature(object = "lmer2"),
	  function(object, REML = 0, ...) {
	      rr <- as(.Call(lmer2_sigma, object, REML)^2 *
                       crossprod(.Call(lmer2_vcov, object)), "dpoMatrix")
              nms <- object@cnames[[".fixed"]]
              dimnames(rr) <- list(nms, nms)
	      rr@factors$correlation <- as(rr, "corMatrix")
	      rr
	  })

## This is modeled a bit after  print.summary.lm :
printMer2 <- function(x, digits = max(3, getOption("digits") - 3),
                      correlation = TRUE, symbolic.cor = FALSE,
                      signif.stars = getOption("show.signif.stars"), ...)
{
    so <- summary(x)
    REML <- so@dims["isREML"]
    llik <- so@logLik
    dev <- so@deviance
    dims <- x@dims
    glz <- so@isG

    cat(so@methTitle, "\n")
    if (!is.null(x@call$formula))
        cat("Formula:", deparse(x@call$formula),"\n")
    if (!is.null(x@call$data))
        cat("   Data:", deparse(x@call$data), "\n")
    if (!is.null(x@call$subset))
        cat(" Subset:",
            deparse(asOneSidedFormula(x@call$subset)[[2]]),"\n")
    if (is(x, "glmer2"))
        cat(" Family: ", so@family$family, "(",
            so@family$link, " link)\n", sep = "")
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

setMethod("summary", signature(object = "lmer2"),
	  function(object, ...) {
	      fcoef <- fixef(object)
	      vcov <- vcov(object)
	      corF <- vcov@factors$correlation
              dims <- object@dims
	      ## DF <- getFixDF(object)
	      coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd) #, DF = DF)
	      REML <- object@dims["isREML"]
	      llik <- logLik(object, REML)
	      dev <- object@deviance

	      glz <- is(object, "glmer")
	      methTitle <-
		  if (glz)
		      paste("Generalized linear mixed model fit using",
			    switch(object@status["glmm"],
                                   "PQL", "Laplace", "AGQ"))
		  else paste("Linear mixed-effects model fit by",
			     if(REML) "REML" else "maximum likelihood")

	      AICframe <- {
		  if (glz)
		      data.frame(AIC = AIC(llik), BIC = BIC(llik),
				 logLik = c(llik),
				 deviance = -2*llik,
				 row.names = "")
		  else
		      data.frame(AIC = AIC(llik), BIC = BIC(llik),
				 logLik = c(llik),
				 MLdeviance = dev["ML"],
				 REMLdeviance = dev["REML"],
				 row.names = "")
	      }
              varcor <- VarCorr(object)
	      REmat <- formatVC(varcor)
              if (is.na(attr(varcor, "sc")))
                  REmat <- REmat[-nrow(REmat), , drop = FALSE]

	      if (nrow(coefs) > 0) {
		  if (dims["famType"] >= 0) {
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

	      new(if(is(object, "glmer")) "summary.glmer" else
                  {if(is(object, "lmer")) "summary.lmer" else "summary.lmer2"},
		  object,
		  isG = glz,
		  methTitle = methTitle,
		  logLik = llik,
		  ngrps = sapply(object@flist, function(x) length(levels(x))),
		  sigma = .Call(lmer2_sigma, object, REML),
		  coefs = coefs,
		  vcov = vcov,
		  REmat = REmat,
		  AICtab= AICframe
		  )
	  })## summary()

## Extract the log-likelihood or restricted log-likelihood
setMethod("logLik", signature(object="lmer2"),
	  function(object, REML = NULL, ...) {
	      dims <- object@dims
              if (is.null(REML) || is.na(REML[1]))
                  REML <- object@dims["isREML"]
	      val <- -deviance(object, REML = REML)/2
	      attr(val, "nall") <- attr(val, "nobs") <- dims["n"]
	      attr(val, "df") <-
                  dims["p"] + length(.Call(lmer2_getPars, object))
	      attr(val, "REML") <-  as.logical(REML)
	      class(val) <- "logLik"
	      val
	  })

## Extract the deviance
setMethod("deviance", signature(object="lmer2"),
	  function(object, REML = NULL, ...) {
              if (missing(REML) || is.null(REML) || is.na(REML[1]))
                  REML <- object@dims["isREML"]
              object@deviance[ifelse(REML, "REML", "ML")]
          })

# Create the VarCorr object of variances and covariances
setMethod("VarCorr", signature(x = "lmer2"),
	  function(x, REML = NULL, ...)
      {
	  sc <- .Call(lmer2_sigma, x, REML)
	  cnames <- x@cnames
	  ans <- x@ST
          for (i in seq(along = ans)) {
              ai <- ans[[i]]
              dm <- dim(ai)
              if (dm[1] < 2) {
                  el <- (sc * ai)^2
              } else {
                  dd <- diag(ai)
                  diag(ai) <- rep(1, dm[1])
                  el <- sc^2 * crossprod(dd * t(ai))
              }
              el <- as(el, "dpoMatrix")
              el@Dimnames <- list(cnames[[i]], cnames[[i]])
	      el@factors$correlation <- as(el, "corMatrix")
	      ans[[i]] <- el
	  }
	  attr(ans, "sc") <- sc
	  ans
      })

# Create the VarCorr object of variances and covariances
setMethod("VarCorr", signature(x = "nlmer"),
	  function(x, ...)
      {
	  sc <- sqrt(sum(x@deviance[c("bqd", "Sdr")])/x@dims["n"])
	  cnames <- x@cnames
	  ans <- x@ST
          for (i in seq(along = ans)) {
              ai <- ans[[i]]
              dm <- dim(ai)
              if (dm[1] < 2) {
                  el <- (sc * ai)^2
              } else {
                  dd <- diag(ai)
                  diag(ai) <- rep(1, dm[1])
                  el <- sc^2 * crossprod(dd * t(ai))
              }
              el <- as(el, "dpoMatrix")
              el@Dimnames <- list(cnames[[i]], cnames[[i]])
	      el@factors$correlation <- as(el, "corMatrix")
	      ans[[i]] <- el
	  }
	  attr(ans, "sc") <- sc
	  ans
      })

setMethod("print", "lmer2", printMer2)
setMethod("show", "lmer2", function(object) printMer2(object))

## Methods for "summary.*" objects:
setMethod("vcov", signature(object = "summary.lmer2"),
	  function(object) object@vcov)
setMethod("logLik", signature(object = "summary.lmer2"),
	  function(object) object@logLik)
setMethod("deviance", signature(object = "summary.lmer2"),
 	  function(object) object@deviance)
setMethod("summary", signature(object = "summary.lmer2"),
          function(object) object)
setMethod("anova", signature(object = "lmer2"),
	  function(object, ...)
      {
	  mCall <- match.call(expand.dots = TRUE)
	  dots <- list(...)
	  modp <- if (length(dots))
	      sapply(dots, is, "lmer2") | sapply(dots, is, "lm") else logical(0)
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
				check.names = FALSE)
	      class(val) <- c("anova", class(val))
	      attr(val, "heading") <-
		  c(header, "Models:",
		    paste(names(mods),
			  unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
			  sep = ": "))
	      return(val)
	  }
	  else { ## ------ single model ---------------------
	      foo <- object
	      ss <- foo@rXy^2
	      ssr <- exp(foo@devComp["logryy2"])
	      names(ss) <- object@cnames[[".fixed"]]
	      asgn <- attr(foo@X, "assign")
	      terms <- foo@terms
	      nmeffects <- attr(terms, "term.labels")
	      if ("(Intercept)" %in% names(ss))
		  nmeffects <- c("(Intercept)", nmeffects)
	      ss <- unlist(lapply(split(ss, asgn), sum))
	      df <- unlist(lapply(split(asgn,  asgn), length))
	      #dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
	      ms <- ss/df
	      #f <- ms/(ssr/dfr)
	      #P <- pf(f, df, dfr, lower.tail = FALSE)
	      #table <- data.frame(df, ss, ms, dfr, f, P)
	      table <- data.frame(df, ss, ms)
	      dimnames(table) <-
		  list(nmeffects,
#			c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
		       c("Df", "Sum Sq", "Mean Sq"))
	      if ("(Intercept)" %in% nmeffects)
		  table <- table[-match("(Intercept)", nmeffects), ]
	      attr(table, "heading") <- "Analysis of Variance Table"
	      class(table) <- c("anova", "data.frame")
	      table
	  }
      })

## Temporary function to convert the ST representation of the
## relative variance-covariance matrix returned by lmer2 into the
## Omega representation required by lmer
ST2Omega <- function(ST)
{
    if (nrow(ST) == 1) return(as(1/ST^2, "dpoMatrix"))
    dd <- diag(ST)
    T <- as(ST, "dtrMatrix")
    T@diag <- "U"
    crossprod(solve(T)/dd)
}

## Extract the random effects
setMethod("ranef", signature(object = "lmer2"),
	  function(object, postVar = FALSE, ...) {
	      ans <- new("ranef.lmer",
                         lapply(.Call(lmer2_ranef, object),
                                data.frame, check.names = FALSE))
              names(ans) <- names(object@flist)
              if (postVar) {
                  pV <- .Call(lmer2_postVar, object)
                  for (i in seq(along = ans))
                      attr(ans[[i]], "postVar") <- pV[[i]]
              }
              ans
	  })

## Generate a Markov chain Monte Carlo sample from the posterior distribution
## of the parameters in a linear mixed model
setMethod("mcmcsamp", signature(object = "lmer2"),
	  function(object, n = 1, verbose = FALSE, saveb = FALSE,
		   trans = TRUE, deviance = FALSE, ...)
      {
          ans <- t(.Call(lmer2_MCMCsamp, object, saveb, n,
                         trans, verbose, deviance))
	  attr(ans, "mcpar") <- as.integer(c(1, n, 1))
	  class(ans) <- "mcmc"
	  mcmccompnames(ans, object, saveb, trans,
			glmer=FALSE, deviance=deviance)
      })

setMethod("simulate", "lmer2",
          function(object, nsim = 1, seed = NULL, ...)
      {
	  if(!is.null(seed)) set.seed(seed)
	  if(!exists(".Random.seed", envir = .GlobalEnv))
	      runif(1)		     # initialize the RNG if necessary
          RNGstate <- .Random.seed
          dims <- object@dims
          ans <- array(0, c(dims["n"], nsim))
          fe <- fixef(object)
          re <- ranef(object)
          nc <- sapply(re, ncol)
          nr <- sapply(re, nrow)
          sz <- nc * nr
          vc <- VarCorr(object)

          cholL <- lapply(vc, chol)
          n <- object@dims["n"]
          for (i in seq_len(nsim))
              ans[, i] <- crossprod(object@ZXyt,
                                    c(unlist(lapply(seq_along(re), function(k)
                                                    (t(cholL[[k]]) %*%
                                                     matrix(rnorm(sz[k]),
                                                            nc = nr[k]))@x)),
                                      fe, 0))@x
          ans + rnorm(prod(dim(ans)), sd = attr(vc, "sc"))
      })

## Remove the intercept row from a transposed model matrix
## first checking that it is not the only column
rmIntr <- function(mm)
{
    if (is.na(irow <- match("(Intercept)", rownames(mm)))) return(mm)
    if (ncol(mm) < 2)
        stop("lhs of a random-effects term cannot be an intercept only")
    mm[-irow, , drop = FALSE]
}


## Fit a nonlinear mixed-effects model
nlmer <- function(formula, data,
                  control = list(), start = NULL, verbose = FALSE,
                  subset, weights, na.action, contrasts = NULL,
                  model = TRUE, ...)
{
    mc <- match.call()
    formula <- as.formula(formula)
    if (length(formula) < 3) stop("formula must be a 3-part formula")
    nlform <- as.formula(formula[[2]])
    if (length(nlform) < 3)
        stop("formula must be a 3-part formula")
    nlmod <- as.call(nlform[[3]])

    cv <- do.call("lmerControl", control)
    if (missing(verbose)) verbose <- cv$msVerbose
    if (is.numeric(start)) start <- list(fixed = start)
    s <- length(pnames <- names(start$fixed))
    stopifnot(length(start$fixed) > 0, s > 0,
              inherits(data, "data.frame"), nrow(data) > 1)
    if (any(pnames %in% names(data)))
        stop("parameter names must be distinct from names of the variables in data")
    anms <- all.vars(nlmod)
    if (!all(pnames %in% anms))
        stop("not all parameter names are used in the nonlinear model expression")

    if (!length(vnms <- setdiff(anms, pnames)))
        stop("there are no variables used in the nonlinear model expression")
    ## create a frame in which to evaluate the factor list
    fr <- lmerFrames(mc,
                     eval(substitute(foo ~ bar,
                                     list(foo = nlform[[2]],
                                          bar = subnms(formula[[3]], lapply(pnames, as.name))))),
                     contrasts, vnms)
    mf <- fr$mf
    attr(mf, "terms") <- NULL
    env <- new.env()
    lapply(names(mf), function(nm) assign(nm, env = env, mf[[nm]]))
    n <- nrow(mf)
    lapply(pnames, function(nm) assign(nm, env = env, rep(start$fixed[[nm]], length.out = n)))

    mf <- mf[rep(seq_len(nrow(data)), s), ]
    row.names(mf) <- seq_len(nrow(mf))
    ss <- rep.int(nrow(data), s)
    for (nm in pnames) mf[[nm]] <- rep.int(as.numeric(nm == pnames), ss)
                                        # factor list and model matrices
    FL <- lmerFactorList(substitute(foo ~ bar, list(foo = nlform[[2]], bar = formula[[3]])),
                         mf, -1)
    FL$Ztl <- lapply(FL$Ztl, rmIntr)    # Remove any intercept columns

    Xt <- t(Matrix(as.matrix(mf[,pnames]), sparse = TRUE))
    xnms <- colnames(fr$X)
    if (!is.na(icol <- match("(Intercept)",xnms))) xnms <- xnms[-icol]
    Xt@Dimnames[[2]] <- NULL
### FIXME: The only times there would be additional columns in the
### fixed effects would be as interactions with parameter names and
### they must be constructed differently
    if (length(xnms) > 0)
        Xt <- rBind(Xt, t(Matrix(fr$X[rep.int(seq_len(n), s), xnms, drop = FALSE])))

    wts <- fr$weights
    if (is.null(wts)) wts <- numeric(0)
    Ztl1 <- lapply(with(FL, .Call(Ztl_sparse, fl, Ztl)), drop0)
    Gp <- unname(c(0L, cumsum(unlist(lapply(Ztl1, nrow)))))
    Zt <- do.call(rBind, Ztl1)
    attr(fr$mf, "terms") <- NULL
    fixef <- unlist(start$fixed)
    storage.mode(fixef) <- "double"
    val <- .Call(nlmer_create, env, nlmod, fr$mf, pnames, call = mc,
                 FL$fl, Xt, Zt, unname(fr$Y), wts,
                 cnames = lapply(FL$Ztl, rownames), Gp = Gp,
                 fixef = fixef)
    .Call(nlmer_optimize, val, verbose)
    .Call(nlmer_update_ranef, val)
    val
}

setMethod("show", "nlmer", function(object)
      {
          dims <- object@dims
          cat("Nonlinear mixed model fit by Laplace\n")
          if (!is.null(object@call$formula))
              cat("Formula:", deparse(object@call$formula),"\n")
          if (!is.null(object@call$data))
              cat("   Data:", deparse(object@call$data), "\n")
          if (!is.null(object@call$subset))
              cat(" Subset:",
                  deparse(asOneSidedFormula(object@call$subset)[[2]]),"\n")

          cat("Random effects:\n")
          print(formatVC(VarCorr(object)), quote = FALSE,
                digits = max(3, getOption("digits") - 3))

          cat(sprintf("Number of obs: %d, groups: ", dims["n"]))
          ngrps <- sapply(object@flist, function(x) length(levels(x)))
          cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
          cat("\n")
          cat("\nFixed effects:\n")
          print(object@fixef)
          invisible(object)
      })


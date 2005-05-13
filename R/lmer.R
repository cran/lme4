## Methods for lmer and for the objects that it produces

## Some utilities

Lind <- function(i,j) {
    if (i < j) stop(paste("Index i=", i,"must be >= index j=", j))
    ((i - 1) * i)/2 + j
}

Dhalf <- function(from) {
    D <- from@D
    nf <- length(D)
    Gp <- from@Gp
    res <- array(0, rep(Gp[nf+1],2))
    for (i in 1:nf) {
        DD <- D[[i]]
        dd <- dim(DD)
        for (k in 1:dd[3]) {
            mm <- array(DD[ , , k], dd[1:2])
            base <- Gp[i] + (k - 1)*dd[1]
            res[cbind(c(base + row(mm)), c(base + col(mm)))] <- c(mm)
        }
    }
    res
}

lmerControl <-                            # Control parameters for lmer
  function(maxIter = 50,
           msMaxIter = 50,
           tolerance = sqrt((.Machine$double.eps)),
           niterEM = 20,
           msTol = sqrt(.Machine$double.eps),
           msVerbose = getOption("verbose"),
           PQLmaxIt = 20,
           EMverbose = getOption("verbose"),
           analyticGradient = TRUE,
           analyticHessian=FALSE)
{
    list(maxIter = maxIter,
         msMaxIter = msMaxIter,
         tolerance = tolerance,
         niterEM = niterEM,
         msTol = msTol,
         msVerbose = msVerbose,
         PQLmaxIt = PQLmaxIt,
         EMverbose=EMverbose,
         analyticHessian=analyticHessian,
         analyticGradient=analyticGradient)
}

setMethod("lmer", signature(formula = "formula", family = "missing"),
          function(formula, data, family,
                   method = c("REML", "ML", "PQL", "Laplace", "AGQ"),
                   control = list(),
                   subset, weights, na.action, offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
                                        # match and check parameters
          method <- match.arg(method)
          if (method %in% c("PQL", "Laplace", "AGQ")) {
              warning(paste('Argument method = "', method,
                            '" is not meaningful for a linear mixed model.\n',
                            'Using method = "REML".\n', sep = ''))
              method <- "REML"
          }
          controlvals <- do.call("lmerControl", control)
          controlvals$REML <- REML <- method == "REML"

          if (length(formula) < 3) stop("formula must be a two-sided formula")

          mf <- match.call()           # create the model frame as frm
          m <- match(c("data", "subset", "weights", "na.action", "offset"),
                     names(mf), 0)
          mf <- mf[c(1, m)]
          mf[[1]] <- as.name("model.frame")
          frame.form <- subbars(formula)
          environment(frame.form) <- environment(formula)
          mf$formula <- frame.form
          mf$drop.unused.levels <- TRUE
          frm <- eval(mf, parent.frame())
          
          ## grouping factors and model matrices for random effects
          bars <- findbars(formula[[3]])
          random <-
              lapply(bars,
                     function(x) list(model.matrix(eval(substitute(~term,
                                                                   list(term=x[[2]]))),
                                                   frm),
                                      eval(substitute(as.factor(fac)[,drop = TRUE],
                                                      list(fac = x[[3]])), frm)))
          names(random) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
          
          ## order factor list by decreasing number of levels
          nlev <- sapply(random, function(x) length(levels(x[[2]])))
          if (any(diff(nlev) > 0)) {
              random <- random[rev(order(nlev))]
          }
          fixed.form <- nobars(formula)
          if (!inherits(fixed.form, "formula")) fixed.form <- ~ 1 # default formula
          Xmat <- model.matrix(fixed.form, frm)
          mmats <- c(lapply(random, "[[", 1),
                     .fixed = list(cbind(Xmat, .response = model.response(frm))))
          obj <- .Call("lmer_create", lapply(random, "[[", 2),
                       mmats, PACKAGE = "Matrix")
          slot(obj, "frame") <- frm
          slot(obj, "terms") <- attr(model.frame(fixed.form, data), "terms")
          slot(obj, "assign") <- attr(Xmat, "assign")
          slot(obj, "call") <- match.call()
          slot(obj, "REML") <- REML
          rm(Xmat)
          .Call("lmer_initial", obj, PACKAGE="Matrix")
          .Call("lmer_ECMEsteps", obj, 
                controlvals$niterEM,
                controlvals$REML,
                controlvals$EMverbose,
                PACKAGE = "Matrix")
          LMEoptimize(obj) <- controlvals
          slot(obj, "residuals") <-
              unname(model.response(frm) -
                     (slot(obj, "fitted") <-
                      .Call("lmer_fitted", obj, mmats, TRUE, PACKAGE = "Matrix")))
          obj
      })

setReplaceMethod("LMEoptimize", signature(x="lmer", value="list"),
                 function(x, value)
             {
                 if (value$msMaxIter < 1) return(x)
                 st <- ccoef(x)         # starting values
                 nc <- x@nc
                 nc <- nc[1:(length(nc) - 2)]
                 constr <- unlist(lapply(nc, function(k) 1:((k*(k+1))/2) <= k))
                 fn <- function(pars) {
                     ccoef(x) <- pars
                     deviance(x, REML = value$REML)
                 }
                 gr <- if (value$analyticGradient)
                     function(pars) {
                         if (!identical(TRUE,all.equal(pars, ccoef(x)))) ccoef(x) <- pars
                         grad <- gradient(x, REML = value$REML, unconst = TRUE)
                         grad[constr] <- -grad[constr]/pars[constr]
                         grad
                     } else NULL
                 optimRes <- optim(st, fn, gr,
                                   method = "L-BFGS-B",
                                   lower = ifelse(constr, 1e-10, -Inf),
                                   control = list(maxit = value$msMaxIter,
                                   trace = as.integer(value$msVerbose)))
                 if (optimRes$convergence != 0) {
                     warning(paste("optim returned message",optimRes$message,"\n"))
                 }
                 ccoef(x) <- optimRes$par
                 return(x)
             })

setMethod("ranef", signature(object = "lmer"),
          function(object, accumulate = FALSE, ...) {
              val <- new("lmer.ranef",
                         lapply(.Call("lmer_ranef", object, PACKAGE = "Matrix"),
                                data.frame, check.names = FALSE),
                         varFac = object@bVar,
                         stdErr = .Call("lmer_sigma", object,
                         object@REML, PACKAGE = "Matrix"))
              if (!accumulate || length(val@varFac) == 1) return(val)
              ## check for nested factors
              L <- object@L
              if (any(sapply(seq(a = val), function(i) length(L[[Lind(i,i)]]@i))))
                  error("Require nested grouping factors to accumulate random effects")
              val
          })

setMethod("fixef", signature(object = "lmer"),
          function(object, ...) {
              val <- .Call("lmer_fixef", object, PACKAGE = "Matrix")
              val[-length(val)]
          })

setMethod("fixef", signature(object = "glmer"),
          function(object, ...) {
              object@fixed
          })

setMethod("VarCorr", signature(x = "lmer"),
          function(x, REML = TRUE, useScale = TRUE, ...) {
              val <- .Call("lmer_variances", x, PACKAGE = "Matrix")
              for (i in seq(along = val)) {
                  dimnames(val[[i]]) = list(x@cnames[[i]], x@cnames[[i]])
                  val[[i]] = as(as(val[[i]], "pdmatrix"), "corrmatrix")
              }
              new("VarCorr",
                  scale = .Call("lmer_sigma", x, REML, PACKAGE = "Matrix"),
                  reSumry = val,
                  useScale = useScale)
          })

setMethod("gradient", signature(x = "lmer"),
          function(x, REML, unconst, ...)
          .Call("lmer_gradient", x, REML, unconst, PACKAGE = "Matrix"))

setMethod("summary", signature(object = "lmer"),
          function(object, ...)
          new("summary.lmer", object, useScale = TRUE,
              showCorrelation = TRUE,
              method = if (object@REML) "REML" else "ML",
              family = gaussian(),
              logLik = logLik(object),
              fixed = fixef(object)))

## FIXME: glmm-s with scale not handled 
setMethod("summary", signature(object = "glmer"),
          function(object, ...)
          new("summary.lmer", object, useScale = FALSE,
              showCorrelation = TRUE,
              method = object@method,
              family = object@family,
              logLik = logLik(object),
              fixed = fixef(object)))

setMethod("show", signature(object = "lmer"),
          function(object)
          show(new("summary.lmer", object, useScale = TRUE,
                   showCorrelation = FALSE,
                   method = if (object@REML) "REML" else "ML",
                   family = gaussian(),
                   logLik = logLik(object),
                   fixed = fixef(object)))
          )

setMethod("show", signature(object = "glmer"),
          function(object)
          show(new("summary.lmer", object, useScale = FALSE,
                   showCorrelation = FALSE,
                   method = object@method,
                   family = object@family,
                   logLik = logLik(object),
                   fixed = fixef(object)))
          )

setMethod("show", "summary.lmer",
          function(object) {
              fcoef <- object@fixed
              useScale <- object@useScale
              corF <- as(as(vcov(object, useScale = useScale), "pdmatrix"),
                         "corrmatrix")
              DF <- getFixDF(object)
              coefs <- cbind(fcoef, corF@stdDev, DF)
              nc <- object@nc
              dimnames(coefs) <-
                  list(names(fcoef), c("Estimate", "Std. Error", "DF"))
                            digits <- max(3, getOption("digits") - 2)
              REML <- length(object@REML) > 0 && object@REML[1]
              llik <- object@logLik
              dev <- object@deviance
              
              rdig <- 5
              if (glz <- !(object@method %in% c("REML", "ML"))) {
                  cat(paste("Generalized linear mixed model fit using",
                            object@method, "\n"))
              } else {
                  cat("Linear mixed-effects model fit by ")
                  cat(if(object@REML) "REML\n" else "maximum likelihood\n")
              }
              if (!is.null(object@call$formula)) {
                  cat("Formula:", deparse(object@call$formula),"\n")
              }
              if (!is.null(object@call$data)) {
                  cat("   Data:", deparse(object@call$data), "\n")
              }
              if (!is.null(object@call$subset)) {
                  cat(" Subset:",
                      deparse(asOneSidedFormula(object@call$subset)[[2]]),"\n")
              }
              if (glz) {
                  cat(" Family:", object@family$family, "(",
                      object@family$link, "link)\n")
                  print(data.frame(AIC = AIC(llik), BIC = BIC(llik),
                               logLik = c(llik),
                               deviance = -2*llik,
                               row.names = ""))
              } else {
                  print(data.frame(AIC = AIC(llik), BIC = BIC(llik),
                               logLik = c(llik),
                               deviance = dev["ML"],
                               REMLdeviance = dev["REML"],
                               row.names = ""))
              }
              cat("Random effects:\n")
              show(VarCorr(object))
              ngrps <- lapply(object@flist, function(x) length(levels(x)))
              cat(sprintf("# of obs: %d, groups: ", object@nc[length(object@nc)]))
              cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
              cat("\n")
              if (!useScale)
                  cat("\nEstimated scale (compare to 1) ",
                      .Call("lmer_sigma", object, object@REML, PACKAGE = "Matrix"),
                      "\n")
              if (nrow(coefs) > 0) {
                  if (useScale) {
                      stat <- coefs[,1]/coefs[,2]
                      pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
                      nms <- colnames(coefs)
                      coefs <- cbind(coefs, stat, pval)
                      colnames(coefs) <- c(nms, "t value", "Pr(>|t|)")
                  } else {
                      coefs <- coefs[, 1:2, drop = FALSE]
                      stat <- coefs[,1]/coefs[,2]
                      pval <- 2*pnorm(abs(stat), lower = FALSE)
                      nms <- colnames(coefs)
                      coefs <- cbind(coefs, stat, pval)
                      colnames(coefs) <- c(nms, "z value", "Pr(>|z|)")
                  }
                  cat("\nFixed effects:\n")
                  printCoefmat(coefs, tst.ind = 4, zap.ind = 3)
                  if (length(object@showCorrelation) > 0 && object@showCorrelation[1]) {
                      rn <- rownames(coefs)
                      dimnames(corF) <- list(
                                               abbreviate(rn, minlen=11),
                                               abbreviate(rn, minlen=6))
                      if (!is.null(corF)) {
                          p <- NCOL(corF)
                          if (p > 1) {
                              cat("\nCorrelation of Fixed Effects:\n")
                              corF <- format(round(corF, 3), nsmall = 3)
                              corF[!lower.tri(corF)] <- ""
                              print(corF[-1, -p, drop=FALSE], quote = FALSE)
                          }
                      }
                  }
              }
              invisible(object)
          })


setMethod("lmer", signature(formula = "formula"),
          function(formula, family, data,
                   method = c("REML", "ML", "PQL", "Laplace", "AGQ"),
                   control = list(),
                   subset, weights, na.action, offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
          if (missing(method)) {
              method <- "PQL"
          } else {
              method <- match.arg(method)
              if (method == "ML") method <- "PQL"
              if (method == "REML") 
                  warning(paste('Argument method = "REML" is not meaningful',
                                'for a generalized linear mixed model.',
                                '\nUsing method = "PQL".\n'))
          }
##           if (method %in% c("Laplace", "AGQ"))
##               stop("Laplace and AGQ methods not yet implemented")
          if (method %in% c("AGQ"))
              stop("AGQ method not yet implemented")
          gVerb <- getOption("verbose")
                                        # match and check parameters
          controlvals <- do.call("lmerControl", control)
          controlvals$REML <- FALSE
          if (length(formula) < 3) stop("formula must be a two-sided formula")
          ## initial glm fit
          mf <- match.call()            
          m <- match(c("family", "data", "subset", "weights",
                       "na.action", "offset"),
                     names(mf), 0)
          mf <- mf[c(1, m)]
          mf[[1]] <- as.name("glm")
          fixed.form <- nobars(formula)
          if (!inherits(fixed.form, "formula")) fixed.form <- ~ 1 # default formula
          environment(fixed.form) <- environment(formula)
          mf$formula <- fixed.form
          mf$x <- mf$model <- mf$y <- TRUE
          glm.fit <- eval(mf, parent.frame())
          family <- glm.fit$family
          ## Note: offset is on the linear scale
          offset <- glm.fit$offset
          if (is.null(offset)) offset <- 0
          weights <- sqrt(abs(glm.fit$prior.weights))
          ## initial 'fitted' values on linear scale
          etaold <- eta <- glm.fit$linear.predictors
          
          ## evaluation of model frame
          mf$x <- mf$model <- mf$y <- mf$family <- NULL
          mf$drop.unused.levels <- TRUE
          this.form <- subbars(formula)
          environment(this.form) <- environment(formula)
          mf$formula <- this.form
          mf[[1]] <- as.name("model.frame")
          frm <- eval(mf, parent.frame())
          
          ## grouping factors and model matrices for random effects
          bars <- findbars(formula[[3]])
          random <-
              lapply(bars,
                     function(x) list(model.matrix(eval(substitute(~term,
                                                                   list(term=x[[2]]))),
                                                   frm),
                                      eval(substitute(as.factor(fac)[,drop = TRUE],
                                                      list(fac = x[[3]])), frm)))
          names(random) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
          
          ## order factor list by decreasing number of levels
          nlev <- sapply(random, function(x) length(levels(x[[2]])))
          if (any(diff(nlev) > 0)) {
              random <- random[rev(order(nlev))]
          }
          mmats <- c(lapply(random, "[[", 1),
                     .fixed = list(cbind(glm.fit$x, .response = glm.fit$y)))
          ## FIXME: Use Xfrm and Xmat to get the terms and assign
          ## slots, pass these to lmer_create, then destroy Xfrm, Xmat, etc.
          obj <- .Call("lmer_create", lapply(random, "[[", 2),
                       mmats, PACKAGE = "Matrix")
          slot(obj, "frame") <- frm
          slot(obj, "terms") <- attr(glm.fit$model, "terms")
          slot(obj, "assign") <- attr(glm.fit$x, "assign")
          slot(obj, "call") <- match.call()
          slot(obj, "REML") <- FALSE
          rm(glm.fit)
          .Call("lmer_initial", obj, PACKAGE="Matrix")
          mmats.unadjusted <- mmats
          mmats[[1]][1,1] <- mmats[[1]][1,1]
          conv <- FALSE
          firstIter <- TRUE
          msMaxIter.orig <- controlvals$msMaxIter
          responseIndex <- ncol(mmats$.fixed)

          for (iter in seq(length = controlvals$PQLmaxIt))
          {
              mu <- family$linkinv(eta)
              dmu.deta <- family$mu.eta(eta)
              ## weights (note: weights is already square-rooted)
              w <- weights * dmu.deta / sqrt(family$variance(mu))
              ## adjusted response (should be comparable to X \beta, not including offset
              z <- eta - offset + (mmats.unadjusted$.fixed[, responseIndex] - mu) / dmu.deta
              .Call("nlme_weight_matrix_list",
                    mmats.unadjusted, w, z, mmats, PACKAGE="Matrix")
              .Call("lmer_update_mm", obj, mmats, PACKAGE="Matrix")
              if (firstIter) {
                  .Call("lmer_initial", obj, PACKAGE="Matrix")
                  if (gVerb) cat(" PQL iterations convergence criterion\n")
              }
              .Call("lmer_ECMEsteps", obj, 
                    controlvals$niterEM,
                    FALSE,
                    controlvals$EMverbose,
                    PACKAGE = "Matrix")
              LMEoptimize(obj) <- controlvals
              eta[] <- offset + ## FIXME: should the offset be here ?
                  .Call("lmer_fitted", obj,
                        mmats.unadjusted, TRUE, PACKAGE = "Matrix")
              crit <- max(abs(eta - etaold)) / (0.1 + max(abs(eta)))
              if (gVerb) cat(sprintf("%03d: %#11g\n", as.integer(iter), crit))
              ## use this to determine convergence
              if (crit < controlvals$tolerance) {
                  conv <- TRUE
                  break
              }
              etaold[] <- eta

              ## Changing number of iterations on second and
              ## subsequent iterations.
              if (firstIter)
              {
                  controlvals$niterEM <- 2
                  controlvals$msMaxIter <- 10
                  firstIter <- FALSE
              }
          }
          if (!conv) warning("IRLS iterations for glmm did not converge")
          controlvals$msMaxIter <- msMaxIter.orig


###          if (TRUE) ## Laplace
###          {
          ## Need to optimize L(theta, beta) using Laplace approximation

          ## Things needed for that:
          ##
          ## 1. reduced ssclme object, offset, weighted model matrices
          ## 2. facs, reduced model matrices

          ## Of these, those in 2 will be fixed given theta and beta,
          ## and can be thought of arguments to the L(theta, beta)
          ## function. However, the ones in 1 will have the same
          ## structure. So the plan is to pre-allocate them and pass
          ## them in too so they can be used without creating/copying
          ## them more than once


          ## reduced ssclme

          reducedObj <- .Call("lmer_collapse", obj, PACKAGE = "Matrix")
          reducedMmats.unadjusted <- mmats.unadjusted
          reducedMmats.unadjusted$.fixed <-
              reducedMmats.unadjusted$.fixed[, responseIndex, drop = FALSE]
          reducedMmats <- mmats
          reducedMmats$.fixed <-
              reducedMmats$.fixed[, responseIndex, drop = FALSE]

          ## define function that calculates bhats given theta and beta 

          bhat <- 
              function(pars = NULL) # 1:(responseIndex-1) - beta, rest - theta
              {
                  if (is.null(pars))
                  {
                      off <- drop(mmats.unadjusted$.fixed %*%
                                  c(fixef(obj), 0)) + offset
                  }
                  else
                  {
                      .Call("lmer_coefGets",
                            reducedObj,
                            as.double(pars[responseIndex:length(pars)]),
                            TRUE,
                            PACKAGE = "Matrix")
                      off <- drop(mmats.unadjusted$.fixed %*%
                                  c(pars[1:(responseIndex-1)], 0)) + offset
                  }
                  niter <- 20
                  conv <- FALSE

                  eta <- offset + 
                      .Call("lmer_fitted",
                            obj, mmats.unadjusted, TRUE,
                            PACKAGE = "Matrix")
                  etaold <- eta
                  
                  for (iter in seq(length = niter))
                  {
                      mu <- family$linkinv(eta)
                      dmu.deta <- family$mu.eta(eta)
                      w <- weights * dmu.deta / sqrt(family$variance(mu))
                      z <- eta - off + (reducedMmats.unadjusted$.fixed[, 1]
                                        - mu) / dmu.deta 
                      .Call("nlme_weight_matrix_list",
                            reducedMmats.unadjusted, w, z, reducedMmats,
                            PACKAGE="Matrix")
                      .Call("lmer_update_mm",
                            reducedObj, reducedMmats,
                            PACKAGE="Matrix")
                      eta[] <- off + 
                          .Call("lmer_fitted", reducedObj, 
                                reducedMmats.unadjusted, TRUE,
                                PACKAGE = "Matrix")
                      ##cat(paste("bhat Criterion:", max(abs(eta - etaold)) /
                      ##          (0.1 + max(abs(eta))), "\n"))
                      ## use this to determine convergence
                      if (max(abs(eta - etaold)) <
                          (0.1 + max(abs(eta))) * controlvals$tolerance)
                      {
                          conv <- TRUE
                          break
                      }
                      etaold[] <- eta
                      
                  }
                  if (!conv) warning("iterations for bhat did not converge")

                  ## bhat doesn't really need to return anything, we
                  ## just want the side-effect of modifying reducedObj
                  ## In particular, we are interested in
                  ## ranef(reducedObj) and reducedObj@bVar (?). But
                  ## the mu-scale response will be useful for log-lik
                  ## calculations later, so return them anyway

                  invisible(family$linkinv(eta)) 
              }

          ## function that calculates log likelihood (the thing that
          ## needs to get evaluated at each Gauss-Hermite location)
          
          ## log scale ? worry about details later, get the pieces in place
          
          ## this is for the Laplace approximation only. GH is more
          ## complicated 

          devLaplace <- function(pars = NULL)
          {
              ## FIXME: This actually returns half the deviance.
              
              ## gets correct values of bhat and bvars. As a side
              ## effect, mu now has fitted values
              mu <- bhat(pars = pars)

              ## GLM family log likelihood (upto constant ?)(log scale)
              ## FIXME: need to adjust for sigma^2 for appropriate models (not trivial!)

              ## Keep everything on (log) likelihood scale
              
              ## log lik from observations given fixed and random effects
              ## get deviance, then multiply by -1/2 (since deviance = -2 log lik)
              ans <- -sum(family$dev.resids(y = mmats.unadjusted$.fixed[, responseIndex],
                                            mu = mu,
                                            wt = weights^2))/2
              
##               if (is.null(getOption("laplaceinR"))) 
##               {
              ans <- ans +
                  .Call("lmer_laplace_devComp", reducedObj,
                        PACKAGE = "Matrix")
##               }
##               else
##               {
##                   ranefs <- .Call("lmer_ranef", reducedObj, PACKAGE = "Matrix")
##                   ## ans <- ans + reducedObj@devComp[2]/2 # log-determinant of Omega

##                   Omega <- reducedObj@Omega
##                   for (i in seq(along = ranefs))
##                   {
##                       ## contribution for random effects (get it working,
##                       ## optimize later) 
##                       ## symmetrize RE variance
##                       Omega[[i]] <- Omega[[i]] + t(Omega[[i]])
##                       diag(Omega[[i]]) <- diag(Omega[[i]]) / 2

##                       ## want log of `const det(Omega) exp(-1/2 b'
##                       ## Omega b )` i.e., const + log det(Omega) - .5
##                       ## * (b' Omega b)

##                       ## FIXME: need to adjust for sigma^2 for appropriate
##                       ## models (easy).  These are all the b'Omega b,
##                       ## summed as they eventually need to be.  Think of
##                       ## this as sum(rowSums((ranefs[[i]] %*% Omega[[i]])
##                       ## * ranefs[[i]]))

##                       ranef.loglik.det <- nrow(ranefs[[i]]) *
##                           determinant(Omega[[i]], logarithm = TRUE)$modulus/2

## #                      print(ranef.loglik.det)

##                       ranef.loglik.re <-
##                           -sum((ranefs[[i]] %*% Omega[[i]]) * ranefs[[i]])/2

## #                      print(ranef.loglik.re)
                      
##                       ranef.loglik <- ranef.loglik.det + ranef.loglik.re

##                       ## Jacobian adjustment
##                       log.jacobian <- 
##                           sum(log(abs(apply(reducedObj@bVar[[i]],
##                                             3,

##                                             ## next line depends on
##                                             ## whether bVars are variances
##                                             ## or Cholesly factors

##                                             ## function(x) sum(diag(x)))
## ## Was this a bug?                          function(x) sum(diag( La.chol( x ) )))
##                                             function(x) prod(diag( La.chol( x ) )))
##                                       )))

## #                      print(log.jacobian)


##                       ## the constant terms from the r.e. and the final
##                       ## Laplacian integral cancel out both being:
##                       ## ranef.loglik.constant <- 0.5 * length(ranefs[[i]]) * log(2 * base::pi)

##                       ans <- ans + ranef.loglik + log.jacobian
##                   }
##               }
              ## ans is (up to some constant) log of the Laplacian
              ## approximation of the likelihood. Return it's negative
              ## to be minimized

              ##              cat("Parameters: ")
              ##              print(pars)

              ##              cat("Value: ")
              ##              print(as.double(-ans))

              -ans 
          }

          if (method == "Laplace")
          {
###Rprof("/tmp/Laplace-profile.out") # trying to figure out if C-ifying bhat is worthwhile
              ## no analytic gradients or hessians
              optimRes <-
                  optim(fn = devLaplace,
                        par =
                        c(fixef(obj),
                          .Call("lmer_coef",
                                obj,
                                TRUE,
                                PACKAGE = "Matrix")),
                        ## WAS: coef(obj, unconst = TRUE)),
                        method = "BFGS", hessian = TRUE,
                        control = list(trace = getOption("verbose"),
                        reltol = controlvals$msTol,
                        maxit = controlvals$msMaxIter))
              if (optimRes$convergence != 0)
                  warning("optim failed to converge")
              optpars <- optimRes$par
              Hessian <- optimRes$hessian
              
              ##fixef(obj) <- optimRes$par[seq(length = responseIndex - 1)]
              if (getOption("verbose")) {
                  cat(paste("optim convergence code",
                            optimRes$convergence, "\n"))
                  cat("Fixed effects:\n")
                  print(fixef(obj))
                  print(optimRes$par[seq(length = responseIndex - 1)])
                  cat("(Unconstrained) variance coefficients:\n")
                  
                  print(
                        .Call("lmer_coef",
                              obj,
                              TRUE,
                              PACKAGE = "Matrix"))

                  ##coef(obj, unconst = TRUE) <-
                  ##    optimRes$par[responseIndex:length(optimRes$par)]
                  ##print(coef(obj, unconst = TRUE))
                  print( optimRes$par[responseIndex:length(optimRes$par)] )
              }
                  
              ## need to calculate likelihood.  also need to store
              ## new estimates of fixed effects somewhere
              ## (probably cannot update standard errors)
###Rprof(NULL)
          }
          else
          {
              optpars <-
                  c(fixef(obj),
                    .Call("lmer_coef",
                          obj,
                          TRUE,
                          PACKAGE = "Matrix"))
              Hessian <- new("matrix")
          }


          ## Before finishing, we need to call devLaplace with the
          ## optimum pars to get the final log likelihood (still need
          ## to make sure it's the actual likelihood and not a
          ## multiple). This would automatically call bhat() and hence
          ## have the 'correct' random effects in reducedObj.

          loglik <- -devLaplace(optpars)
          ##print(loglik)
          ff <- optpars[1:(responseIndex-1)]
          names(ff) <- names(fixef(obj))

          if (!x) mmats <- list()

###      }

          ##obj

          new("glmer", obj,
              family = family,
              glmmll = loglik,
              method = method,
              fixed = ff)
      })


## calculates degrees of freedom for fixed effects Wald tests
## This is a placeholder.  The answers are generally wrong.  It will
## be very tricky to decide what a 'right' answer should be with
## crossed random effects.

setMethod("getFixDF", signature(object="lmer"),
          function(object, ...)
      {
          nc <- object@nc[-seq(along = object@Omega)]
          p <- nc[1] - 1
          n <- nc[2]
          rep(n-p, p)
      })

setMethod("logLik", signature(object="lmer"),
          function(object, REML = object@REML, ...) {
              val <- -deviance(object, REML = REML)/2
              nc <- object@nc[-seq(a = object@Omega)]
              attr(val, "nall") <- attr(val, "nobs") <- nc[2]
              attr(val, "df") <- nc[1] + length(ccoef(object))
              attr(val, "REML") <- REML 
              class(val) <- "logLik"
              val
          })

setMethod("logLik", signature(object="glmer"),
          function(object, ...) {
              val <- object@glmmll
              nc <- object@nc[-seq(a = object@Omega)]
              attr(val, "nall") <- attr(val, "nobs") <- nc[2]
              attr(val, "df") <- nc[1] + length(ccoef(object))
              ## attr(val, "REML") <- REML 
              class(val) <- "logLik"
              val
          })

setMethod("anova", signature(object = "lmer"),
          function(object, ...)
      {
          mCall <- match.call(expand.dots = TRUE)
          dots <- list(...)
          modp <- logical(0)
          if (length(dots))
              modp <- sapply(dots, inherits, "lmer") | sapply(dots, inherits, "lm")
          if (any(modp)) {              # multiple models - form table
              opts <- dots[!modp]
              mods <- c(list(object), dots[modp])
              names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)], as.character)
              mods <- mods[order(sapply(lapply(mods, logLik, REML = FALSE), attr, "df"))]
              calls <- lapply(mods, slot, "call")
              data <- lapply(calls, "[[", "data")
              if (any(data != data[[1]])) stop("all models must be fit to the same data object")
              header <- paste("Data:", data[[1]])
              subset <- lapply(calls, "[[", "subset")
              if (any(subset != subset[[1]])) stop("all models must use the same subset")
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
          } else {
              foo <- object
              foo@status["factored"] <- FALSE
              .Call("lmer_factor", foo, PACKAGE="Matrix")
              dfr <- getFixDF(foo)
              rcol <- ncol(foo@RXX)
              ss <- foo@RXX[ , rcol]^2
              ssr <- ss[[rcol]]
              ss <- ss[seq(along = dfr)]
              names(ss) <- object@cnames[[".fixed"]][seq(along = dfr)]
              asgn <- foo@assign
              terms <- foo@terms
              nmeffects <- attr(terms, "term.labels")
              if ("(Intercept)" %in% names(ss))
                  nmeffects <- c("(Intercept)", nmeffects)
              ss <- unlist(lapply(split(ss, asgn), sum))
              df <- unlist(lapply(split(asgn,  asgn), length))
              dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
              ms <- ss/df
              f <- ms/(ssr/dfr)
              P <- pf(f, df, dfr, lower.tail = FALSE)
              table <- data.frame(df, ss, ms, dfr, f, P)
              dimnames(table) <-
                  list(nmeffects,
                       c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
              if ("(Intercept)" %in% nmeffects) table <- table[-1,]
              attr(table, "heading") <- "Analysis of Variance Table"
              class(table) <- c("anova", "data.frame")
              table
          }
      })

setMethod("update", signature(object = "lmer"),
          function(object, formula., ..., evaluate = TRUE)
      {
          call <- object@call
          if (is.null(call))
              stop("need an object with call component")
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


setMethod("confint", signature(object = "lmer"),
          function (object, parm, level = 0.95, ...) 
      {
          cf <- fixef(object)
          pnames <- names(cf)
          if (missing(parm)) 
              parm <- seq(along = pnames)
          else if (is.character(parm)) 
              parm <- match(parm, pnames, nomatch = 0)
          a <- (1 - level)/2
          a <- c(a, 1 - a)
          pct <- paste(round(100 * a, 1), "%")
          ci <- array(NA, dim = c(length(parm), 2),
                      dimnames = list(pnames[parm], pct))
          ses <- sqrt(diag(vcov(object)))[parm]
          ci[] <- cf[parm] + ses * t(outer(a, getFixDF(object)[parm], qt))
          ci
      })

setMethod("param", signature(object = "lmer"),
          function(object, unconst = FALSE, ...) {
              .Call("lmer_coef", object, unconst, PACKAGE = "Matrix")
          })

setMethod("deviance", "lmer",
          function(object, REML = NULL, ...) {
              .Call("lmer_factor", object, PACKAGE = "Matrix")
              if (is.null(REML))
                  REML <- if (length(oR <- object@REML)) oR else FALSE
              object@deviance[[ifelse(REML, "REML", "ML")]]
          })


setMethod("deviance", "glmer",
          function(object, ...) {
              -2 * object@glmmll
          })

setMethod("chol", signature(x = "lmer"),
          function(x, pivot = FALSE, LINPACK = pivot) {
              x@status["factored"] <- FALSE # force a decomposition
              .Call("lmer_factor", x, PACKAGE = "Matrix")
          })

setMethod("solve", signature(a = "lmer", b = "missing"),
          function(a, b, ...)
          .Call("lmer_invert", a, PACKAGE = "Matrix")
          )

setMethod("formula", "lmer", function(x, ...) x@call$formula)

setMethod("vcov", signature(object = "lmer"),
          function(object, REML = object@REML, useScale = TRUE,...) {
              sc <- .Call("lmer_sigma", object, REML, PACKAGE = "Matrix")
              rr <- object@RXX
              nms <- object@cnames[[".fixed"]]
              dimnames(rr) <- list(nms, nms)
              nr <- nrow(rr)
              rr <- rr[-nr, -nr, drop = FALSE]
              rr <- rr %*% t(rr)
              if (useScale) {
                  rr = sc^2 * rr
              }
              rr
          })

## Extract the L matrix 
setAs("lmer", "dtTMatrix",
      function(from)
  {
      ## force a refactorization if the factors have been inverted
      if (from@status["inverted"]) from@status["factored"] <- FALSE
      .Call("lmer_factor", from, PACKAGE = "Matrix")
      L <- lapply(from@L, as, "dgTMatrix")
      nf <- length(from@D)
      Gp <- from@Gp
      nL <- Gp[nf + 1]
      Li <- integer(0)
      Lj <- integer(0)
      Lx <- double(0)
      for (i in 1:nf) {
          for (j in 1:i) {
              Lij <- L[[Lind(i, j)]]
              Li <- c(Li, Lij@i + Gp[i])
              Lj <- c(Lj, Lij@j + Gp[j])
              Lx <- c(Lx, Lij@x)
          }
      }
      new("dtTMatrix", Dim = as.integer(c(nL, nL)), i = Li, j = Lj, x = Lx,
          uplo = "L", diag = "U")
  })

## Extract the ZZX matrix
setAs("lmer", "dsTMatrix",
      function(from)
  {
      .Call("lmer_inflate", from, PACKAGE = "Matrix")
      ZZpO <- lapply(from@ZZpO, as, "dgTMatrix")
      ZZ <- lapply(from@ZtZ, as, "dgTMatrix")
      nf <- length(ZZpO)
      Gp <- from@Gp
      nZ <- Gp[nf + 1]
      Zi <- integer(0)
      Zj <- integer(0)
      Zx <- double(0)
      for (i in 1:nf) {
          ZZpOi <- ZZpO[[i]]
          Zi <- c(Zi, ZZpOi@i + Gp[i])
          Zj <- c(Zj, ZZpOi@j + Gp[i])
          Zx <- c(Zx, ZZpOi@x)
          if (i > 1) {
              for (j in 1:(i-1)) {
                  ZZij <- ZZ[[Lind(i, j)]]
                  ## off-diagonal blocks are transposed
                  Zi <- c(Zi, ZZij@j + Gp[j])
                  Zj <- c(Zj, ZZij@i + Gp[i])
                  Zx <- c(Zx, ZZij@x)
              }
          }
      }
      new("dsTMatrix", Dim = as.integer(c(nZ, nZ)), i = Zi, j = Zj, x = Zx,
          uplo = "U")
  })

setMethod("fitted", signature(object = "lmer"),
          function(object, ...)
          napredict(attr(object@frame, "na.action"), object@fitted))

setMethod("residuals", signature(object = "lmer"),
          function(object, ...)
          naresid(attr(object@frame, "na.action"), object@residuals))

setMethod("resid", signature(object = "lmer"),
          function(object, ...) do.call("residuals", c(list(object), list(...))))

setMethod("coef", signature(object = "lmer"),
          function(object, ...)
      {
          fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
          ref <- as(ranef(object), "list")
          names(ref) <- names(object@flist)
          val <- lapply(ref, function(x) fef[rep(1, nrow(x)),])
          for (i in seq(a = val)) {
              refi <- ref[[i]]
              row.names(val[[i]]) <- row.names(refi)
              if (!all(names(refi) %in% names(fef)))
                  stop("unable to align random and fixed effects")
              val[[i]][ , names(refi)] <- val[[i]][ , names(refi)] + refi
          }
          new("lmer.coef", val)
      })

setMethod("plot", signature(x = "lmer.coef"),
          function(x, y, ...)
      {
          varying <- unique(do.call("c",
                                    lapply(x, function(el)
                                           names(el)[sapply(el,
                                                            function(col)
                                                            any(col != col[1]))])))
          gf <- do.call("rbind", lapply(x, "[", j = varying))
          gf$.grp <- factor(rep(names(x), sapply(x, nrow)))
          switch(min(length(varying), 3),
                 qqmath(eval(substitute(~ x | .grp,
                                        list(x = as.name(varying[1])))), gf, ...),
                 xyplot(eval(substitute(y ~ x | .grp,
                                        list(y = as.name(varying[1]),
                                             x = as.name(varying[2])))), gf, ...),
                 splom(~ gf | .grp, ...))
      })

setMethod("plot", signature(x = "lmer.ranef"),
          function(x, y, ...)
      {
          lapply(x, function(x) {
              cn <- lapply(colnames(x), as.name)
              switch(min(ncol(x), 3),
                     qqmath(eval(substitute(~ x, list(x = cn[[1]]))), x, ...),
                     xyplot(eval(substitute(y ~ x, list(y = cn[[1]], x = cn[[2]]))),
                            x, ...),
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

setMethod("GLMM", signature(formula = "missing"),
          function(formula, family, data, random,
                   method = c("PQL", "Laplace"),
                   control = list(),
                   subset,
                   weights,
                   na.action,
                   offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
          nCall = mCall = match.call()
          resp = getResponseFormula(data)[[2]]
          cov = getCovariateFormula(data)[[2]]
          nCall$formula = eval(substitute(resp ~ cov))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "Matrix")
      })

setMethod("GLMM", signature(formula = "formula",
                            data = "groupedData", random = "missing"),
          function(formula, family, data, random,
                   method = c("PQL", "Laplace"),
                   control = list(),
                   subset,
                   weights,
                   na.action,
                   offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
          nCall = mCall = match.call()
          cov = formula[[3]]
          grps = getGroupsFormula(data)[[2]]
          nCall$random = eval(substitute(~ cov | grps))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "Matrix")
      })

setMethod("GLMM", signature(random = "formula"),
          function(formula, family, data, random,
                   method = c("PQL", "Laplace"),
                   control = list(),
                   subset,
                   weights,
                   na.action,
                   offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
          nCall = mCall = match.call()
          cov = getCovariateFormula(random)
          nCall$random <- lapply(getGroupsFormula(random, asList = TRUE),
                                 function(f) cov)
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "Matrix")
      })

setMethod("GLMM", signature(formula = "formula",
                            data = "groupedData",
                            random = "list"),
          function(formula, family, data, random,
                   method = c("PQL", "Laplace"),
                   control = list(),
                   subset,
                   weights,
                   na.action,
                   offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
          nCall = mCall = match.call()
          nCall$data <- data@data
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "Matrix")
      })

make.glm.call <- 
    function (mf, frm) 
{
    m <- match(c("formula", "family", "data", "weights",
                 "subset", "na.action", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("glm")
    ## environment(frm) = environment(formula) ???
    mf$formula = frm
    mf$model = FALSE
    mf$x = FALSE
    mf$y = TRUE
    mf
}

setMethod("GLMM",
          signature(formula = "formula",
                    random = "list"),
          function(formula, family, data, random,
                   method = c("PQL", "Laplace"),
                   control = list(),
                   subset,
                   weights,
                   na.action,
                   offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
          gVerb <- getOption("verbose")
          random <-
              lapply(random,
                     get("formula", pos = parent.frame(), mode = "function"))
          controlvals <- do.call("lmeControl", control)
          controlvals$REML <- FALSE
          method <- match.arg(method)

          ## BEGIN glm fit without random effects

          ## several arguments are handled at this point.
          ## what about

          ##   subset, na.action (arguments to model.frame) 

          ## Are these  already handled ?
          ## The rest could be part of ..., why have them explicitly ?
          ## Do we want possibility of supplying control in glm() ?

          glm.fit <- eval(make.glm.call(match.call(expand.dots = TRUE),
                                        formula), parent.frame())
          family <- glm.fit$family
          offset <- if (is.null(glm.fit$offset)) 0 else glm.fit$offset
          weights <- sqrt(abs(glm.fit$prior.weights))
          ## initial 'fitted' values on linear scale
          eta <- glm.fit$linear.predictors
          ## FIXME: does this include offsets (make sure it does)
          etaold <- eta

          ## END using glm fit results
          ## Note: offset is on the linear scale

          ## FIXME: Not clear how (user specified) offset works. It
          ## doesn't make sense for it to be an offset for the
          ## response on the mu scale, so I'm assuming it's on the
          ## linear predictor scale


          data <- eval(make.mf.call(match.call(expand.dots = FALSE),
                                    formula, random), parent.frame())
          facs <- lapply(names(random),
                         function(x) as.factor(eval(as.name(x), envir = data)))
          names(facs) <- names(random)
          ## order factor list by decreasing number of levels
          ford <- rev(order(sapply(facs, function(fac) length(levels(fac)))))
          if (any(ford != seq(a = ford))) { # re-order both facs and random
              facs <- facs[ford]
              random <- random[ford]
          }
          ## creates model matrices
          mmats.unadjusted <-
              c(lapply(random,
                       function(x) model.matrix(formula(x), data = data)),
                list(.Xy =
                     cbind(model.matrix(formula, data = data),
                           .response = glm.fit$y))) #  WAS: model.response(data)
          responseIndex <- ncol(mmats.unadjusted$.Xy)
          obj <-  ## creates ssclme structure
              .Call("ssclme_create", facs,
                    unlist(lapply(mmats.unadjusted, ncol)),
                    PACKAGE = "Matrix")
          facs = facshuffle(obj, facs)
          obj = obj[[1]]
          mmats <- mmats.unadjusted
          ## the next line is to force a copy of mmats, because we are
          ## going to use both mmats and mmats.unadjusted as arguments
          ## in a .Call where one of them will be modified (don't want
          ## the other to be modified as well)
          mmats[[1]][1,1] <- mmats[[1]][1,1]
          conv <- FALSE
          firstIter <- TRUE
          msMaxIter.orig <- controlvals$msMaxIter

          for (iter in seq(length = controlvals$PQLmaxIt))
          {
              mu <- family$linkinv(eta)
              dmu.deta <- family$mu.eta(eta)
              ## weights (note: weights is already square-rooted)
              w <- weights * dmu.deta / sqrt(family$variance(mu))
              ## adjusted response (should be comparable to X \beta, not including offset
              z <- eta - offset + (mmats.unadjusted$.Xy[, responseIndex] - mu) / dmu.deta
              .Call("nlme_weight_matrix_list",
                    mmats.unadjusted, w, z, mmats, PACKAGE="Matrix")
              .Call("ssclme_update_mm", obj, facs, mmats, PACKAGE="Matrix")
              if (firstIter) {
                  .Call("ssclme_initial", obj, PACKAGE="Matrix")
                  if (gVerb) cat(" PQL iterations convergence criterion\n")
              }
              .Call("ssclme_EMsteps", obj,
                    controlvals$niterEM,
                    FALSE, #controlvals$REML,
                    controlvals$EMverbose,
                    PACKAGE = "Matrix")
              LMEoptimize(obj) = controlvals
              eta[] <- offset + ## FIXME: should the offset be here ?
                  .Call("ssclme_fitted", obj, facs,
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

          reducedObj <- .Call("ssclme_collapse", obj, PACKAGE = "Matrix")
          reducedMmats.unadjusted <- mmats.unadjusted
          reducedMmats.unadjusted$.Xy <-
              reducedMmats.unadjusted$.Xy[, responseIndex, drop = FALSE]
          reducedMmats <- mmats
          reducedMmats$.Xy <-
              reducedMmats$.Xy[, responseIndex, drop = FALSE]

          ## define function that calculates bhats given theta and beta 

          bhat <- 
              function(pars = NULL) # 1:(responseIndex-1) - beta, rest - theta
              {
                  if (is.null(pars))
                  {
                      off <- drop(mmats.unadjusted$.Xy %*%
                                  c(fixef(obj), 0)) + offset
                  }
                  else
                  {
                      .Call("ssclme_coefGets",
                            reducedObj,
                            as.double(pars[responseIndex:length(pars)]),
                            TRUE,
                            PACKAGE = "Matrix")
                      off <- drop(mmats.unadjusted$.Xy %*%
                                  c(pars[1:(responseIndex-1)], 0)) + offset
                  }

                  niter <- 20
                  conv <- FALSE

                  eta <- offset + 
                      .Call("ssclme_fitted", obj, facs,
                            mmats.unadjusted, TRUE, PACKAGE = "Matrix")
                  etaold <- eta
                  
                  for (iter in seq(length = niter))
                  {
                      mu <- family$linkinv(eta)
                      dmu.deta <- family$mu.eta(eta)
                      w <- weights * dmu.deta / sqrt(family$variance(mu))
                      z <- eta - off + (reducedMmats.unadjusted$.Xy[, 1]
                                        - mu) / dmu.deta 
                      .Call("nlme_weight_matrix_list",
                            reducedMmats.unadjusted, w, z, reducedMmats,
                            PACKAGE="Matrix")
                      .Call("ssclme_update_mm", reducedObj, facs, reducedMmats,
                            PACKAGE="Matrix")
                      eta[] <- off + 
                          .Call("ssclme_fitted", reducedObj, facs,
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
              ans <- -sum(family$dev.resids(y = mmats.unadjusted$.Xy[,
                                            responseIndex],
                                            mu = mu,
                                            wt = weights^2))/2
                                                
              ranefs <- ranef(reducedObj)
              # ans <- ans + reducedObj@devComp[2]/2 # log-determinant of Omega
              Omega <- reducedObj@Omega
              for (i in seq(along = ranefs))
              {
                  ## contribution for random effects (get it working,
                  ## optimize later) 
                  ## symmetrize RE variance
                  Omega[[i]] <- Omega[[i]] + t(Omega[[i]])
                  diag(Omega[[i]]) <- diag(Omega[[i]]) / 2

                  ## want log of `const det(Omega) exp(-1/2 b' Omega b )`
                  ## i.e., const + log det(Omega) - .5 * (b' Omega b)
                  ## FIXME: need to adjust for sigma^2 for appropriate models (easy)
                  ## these are all the b'Omega b, summed as they eventually need to be
                  ## think of this as sum(rowSums((ranefs[[i]] %*% Omega[[i]]) * ranefs[[i]]))
                  
                  ranef.loglik.det <- nrow(ranefs[[i]]) *
                      determinant(Omega[[i]], logarithm = TRUE)$modulus/2
                  ranef.loglik.re <- -sum((ranefs[[i]] %*% Omega[[i]]) *
                                          ranefs[[i]])/2
                  ranef.loglik <- ranef.loglik.det + ranef.loglik.re

                  ## Jacobian adjustment
                  log.jacobian <- sum(log(abs(apply(reducedObj@bVar[[i]],
                                                    3,
                                                    function(x) sum(diag(x)))
                                              )))

                  ## the constant terms from the r.e. and the final
                  ## Laplacian integral cancel out both being:
                  ## ranef.loglik.constant <- 0.5 * length(ranefs[[i]]) * log(2 * base::pi)

                  ans <- ans + ranef.loglik + log.jacobian
              }
              ## ans is (up to some constant) log of the Laplacian
              ## approximation of the likelihood. Return it's negative
              ## to be minimized

#              cat("Parameters: ")
#              print(pars)

#              cat("Value: ")
#              print(as.double(-ans))

              -ans 
          }

          if (method == "Laplace")
          {
              cat(paste("Using optimizer", controlvals$optim), "\n")

              ## no analytic gradients or hessians
              if (controlvals$optimizer == "optim")
              {
                  optimRes =
                      optim(fn = devLaplace,
                            par = c(fixef(obj), coef(obj, unconst = TRUE)),
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
                      print(coef(obj, unconst = TRUE))
                      coef(obj, unconst = TRUE) <-
                          optimRes$par[responseIndex:length(optimRes$par)]
                      print(coef(obj, unconst = TRUE))
                  }
              }
              else if (controlvals$optimizer == "nlm")
              {
                  ## not sure what the next few lines are for. Copied from Saikat's code
                  ## typsize <- 1/controlvals$msScale(coef(obj))
                  typsize <- rep(1.0, length(coef(obj, unconst = TRUE)) +
                                 responseIndex - 1)
                  if (is.null(controlvals$nlmStepMax))
                      controlvals$nlmStepMax <-
                          max(100 * sqrt(sum((c(fixef(obj),
                                                coef(obj, unconst = TRUE))/
                                              typsize)^2)), 100)
                  nlmRes =
                      nlm(f = devLaplace, 
                          p = c(fixef(obj), coef(obj, unconst = TRUE)),
                          hessian = TRUE,
                          print.level = if (getOption("verbose")) 2 else 0,
                          steptol = controlvals$msTol,
                          gradtol = controlvals$msTol,
                          stepmax = controlvals$nlmStepMax,
                          typsize=typsize,
                          ## fscale=xval,
                          iterlim = controlvals$msMaxIter)
                  if (nlmRes$code > 3)
                      warning("nlm probably failed to converge")
                  optpars <- nlmRes$estimate
                  Hessian <- nlmRes$hessian
                  
                  if (getOption("verbose")) {
                      cat(paste("nlm convergence code", nlmRes$code, "\n"))
                      cat("Fixed effects:\n")
                      print(fixef(obj))
                      print(nlmRes$estimate[seq(length = responseIndex - 1)])
                      cat("(Unconstrained) variance coefficients:\n")
                      print(coef(obj, unconst = TRUE))
                      coef(obj, unconst = TRUE) <-
                          nlmRes$estimate[responseIndex:
                                          length(nlmRes$estimate)]
                      print(coef(obj, unconst = TRUE))
                  }
              }
              
              ## need to calculate likelihood also need to store new
              ## estimates of fixed effects somewhere (probably cannot
              ## update standard errors)
          } else {
              optpars <- c(fixef(obj), coef(obj, unconst = TRUE))
              Hessian <- new("matrix")
          }


          ## Before finishing, we need to call devLaplace with the
          ## optimum pars to get the final log likelihood (still need
          ## to make sure it's the actual likelihood and not a
          ## multiple). This would automatically call bhat() and hence
          ## have the 'correct' random effects in reducedObj.

          loglik <- devLaplace(optpars)
          ff <- optpars[1:(responseIndex-1)]
          names(ff) <- names(fixef(obj))

          if (!x) mmats <- list()
          ## Things to include in returned object: new ranef
          ## estimates, new parameter estimates (fixed effects and
          ## coefs) (with standard errors from hessian ?) and
          ## (Laplace) approximate log likelihood.

          new("GLMM",
              family = family,
              logLik = -loglik,
              fixef = ff,
              Hessian = Hessian,
              method = method,
              call = match.call(),
              facs = facs,
              x = mmats,
              model = if(model) data else data.frame(list()),
              REML = FALSE,
              rep = if(method == 'Laplace') reducedObj else obj,
              fitted = eta)
      })

setMethod("logLik", signature(object = "GLMM"),
          function(object, ...) {
              val = object@logLik
              rr = object@rep
              nc = rr@nc[-seq(a = rr@Omega)]
              attr(val, "nall") = attr(val, "nobs") = nc[2]
              attr(val, "df") = nc[1] + length(coef(rr))
              attr(val, "REML") = FALSE
              class(val) <- "logLik"
              val
          })

setMethod("deviance", signature(object = "GLMM"),
          function(object, ...) -2 * object@logLik)

setMethod("fixef", signature(object = "GLMM"),
          function(object) object@fixef)

setMethod("ranef", signature(object = "GLMM"),
          function(object, ...)
      {
          object = object@rep
          callGeneric()
      })

setMethod("coef", signature(object = "GLMM"),
          function(object, ...)
      {
          object = object@rep
          callGeneric()
      })

setMethod("show", signature(object = "summary.GLMM"),
          function(object)
      {
          rdig <- 5
          cat("Generalized Linear Mixed Model\n\n")
          cat("Family:", object@family$family, "family with",
              object@family$link, "link\n")
          if (!is.null(object@call$formula)) {
              cat("Fixed:", deparse(object@call$formula),"\n")
          }
          if (!is.null(object@call$data)) {
              cat("Data:", deparse(object@call$data), "\n")
          }
          if (!is.null(object@call$subset)) {
              cat(" Subset:",
                  deparse(asOneSidedFormula(object@call$subset)[[2]]),"\n")
          }
          llik = object@logLik
          print(data.frame(AIC = AIC(llik), BIC = BIC(llik),
                           logLik = c(object@logLik), row.names = ""))
          cat("\n")
          object@re@useScale = FALSE
          object@re@showCorrelation = TRUE
          show(object@re)
          invisible(object)
      })

setMethod("show", signature(object = "GLMM"),
          function(object)
      {
          sumry = summary(object)
          rdig <- 5
          cat("Generalized Linear Mixed Model\n\n")
          cat("Family:", object@family$family, "family with",
              object@family$link, "link\n")
          if (!is.null(object@call$formula)) {
              cat("Fixed:", deparse(object@call$formula),"\n")
          }
          if (!is.null(object@call$data)) {
              cat("Data:", deparse(object@call$data ), "\n")
          }
          if (!is.null(object@call$subset)) {
              cat(" Subset:",
                  deparse(asOneSidedFormula(object@call$subset)[[2]]),"\n")
          }
          cat(paste(" log-likelihood: ", sep = ''), logLik(object), "\n")
          sumry@re@useScale = FALSE
          sumry@re@showCorrelation = FALSE
          saveopt = options(show.signif.stars=FALSE)
          on.exit(options(saveopt))
          show(sumry@re)
          if (object@method != "PQL") { # fix up the fixed effects
              print(t(object@fixef))
          }
          invisible(object)
      })

setMethod("summary", signature(object="GLMM"),
          function(object, ...) {
              llik <- logLik(object)
              resd <- residuals(object, type="pearson")
              if (length(resd) > 5) {
                  resd <- quantile(resd)
                  names(resd) <- c("Min","Q1","Med","Q3","Max")
              }
              re = summary(object@rep, REML = FALSE, useScale = FALSE)
              if (object@method != 'PQL') {
                  hess <- object@Hessian
                  corFixed <- as(as(solve(hess[-nrow(hess),-ncol(hess)]),
                                  "pdmatrix"), "corrmatrix")
                  ## FIXME: change the name of Hessian to info, to make it
                  ## explicit that this is the information matrix
                  fcoef <- object@fixef
                  re@coefficients <- array(c(fcoef, corFixed@stdDev),
                                           c(length(fcoef), 2),
                                           list(names(fcoef),
                                                c("Estimate", "Std. Error")))
                  dimnames(corFixed) <- list(names(fcoef), names(fcoef))
                  re@corFixed <- corFixed
              }
              new("summary.GLMM",
                  family = object@family,
                  call = object@call,
                  logLik = llik,
                  re = re,
                  residuals = resd)
          })

setMethod("glmmStruct", signature(formula = "formula",
                                  random = "list"),
          function(formula, random, family, data, nextraCols, method, ...)
      {

          addToFamily <- function(fam)
          {
              if (method == "Laplace")
                  fam$mu2.eta2 <-
                      if (fam$family %in% c("binomial",
                                            "quasibinomial") &&
                          fam$link == "logit") {
                          function(eta, mu, mu.eta) mu.eta*(1-2*mu)
                      } else if (fam$family %in% c("poisson",
                                                   "quasipoisson") &&
                                 fam$link == "log") {
                          function(eta, mu, mu.eta) mu.eta
                      } else {
                          stop("can not handle link ",
                               fam$link , " for family ",
                               fam$family)
                      }
              fam
          }


          mcall1 <- mcall2 <- match.call()

          mcall1$random <- mcall1$nextraCols <- mcall1$method <- NULL
          mcall1[[1]] <- as.name("glm")
          glmFit <- eval(mcall1, parent.frame())
          rm(mcall1)

          nm <- names(mcall2)[-1]
          keep <-
              is.element(nm, c("data", "subset", "na.action",
                               "xlev", "offset", "weights"))
          for(i in nm[!keep]) mcall2[[i]] <- NULL
          allvars <- c(unlist(lapply(random, all.vars)), names(random))
          mcall2$formula <-
              as.formula(paste(paste(deparse(formula), collapse = ''),
                               paste(allvars, collapse="+"), sep = "+"))
          environment(mcall2$formula) <- environment(formula)
          mcall2$drop.unused.levels <- TRUE
          mcall2[[1]] <- as.name("model.frame")
          data <-  eval(mcall2, parent.frame())
          if (is.null(glmFit$offset))
              glmFit$offset <- 0.0
          data[, attr(attr(data, 'terms'), 'response')] <-
              glmFit$linear.predictor + glmFit$residuals - glmFit$offset

          ans <- new("glmm",
                     reStruct = reStruct(fixed=formula,
                                         random=random,
                                         data=data,
                                         weights=sqrt(abs(glmFit$weights)),
                                         REML=FALSE,
                                         nextraCols=nextraCols),
                     frame = data, call = match.call(),
                     na.action = attr(data, "na.action"),
                     family = addToFamily(glmFit$family),
                     method = method)
          rm(data)
          origOrder <- ans@reStruct@origOrder

          ans@origy <- as.numeric(glmFit$y[origOrder])
          ans@n <- local({
              y <- model.response(glmFit$model, "numeric")
              mt <- attr(glmFit$model, "terms")
              x <- if (!is.empty.model(mt))
                  model.matrix(mt, glmFit$model, glmFit$contrasts)
              x <- as.matrix(x)
              nobs <- NROW(y)
              weights <- glmFit$prior.weights
              ## calculates mustart and may change y and weights and set n (!)
              eval(ans@family$initialize)
              n
          })[origOrder]
          if (length(glmFit$offset) == 1)
              ans@reStruct@offset <- glmFit$offset
          else ans@reStruct@offset <- glmFit$offset[origOrder]
          if (length(glmFit$prior.weights) == 1)
              ans@prior.weights <- glmFit$prior.weights
          else ans@prior.weights <- glmFit$prior.weights[origOrder]
          if (length(glmFit$weights) == 1)
              ans@init.weights <- glmFit$weights
          else ans@init.weights <- glmFit$weights[origOrder]
          ans@init.y <- getResponse(ans@reStruct)
          ans@reStruct@logLik <- as.numeric(NA)
          ans
      })

glmmLa2RespWt <-
    function(fam, eta, origy, w, off)
{
    mu <- fam$linkinv(eta)
    mu.eta.val <- fam$mu.eta(eta)
    list(eta + (origy - mu)/mu.eta.val  - off,
         sqrt(abs(w * mu.eta.val)))
}

glmmLa2Wt2 <-
    function(fam, eta, w)
{
    abs(w)*fam$mu2.eta2(eta, fam$linkinv(eta), fam$mu.eta(eta))
}

glmmLa2LogLikComp <-
    function(x, eta)
{
    fam <- x@family
    origy <- x@origy
    w <- x@prior.weights

    mu <- fam$linkinv(eta)
    nobs <- length(origy)
    aic <- fam$aic(origy, x@n, mu, w, sum(fam$dev.resids(origy, mu, w)))
    ## allow for estimated dispersion
    if (fam$family %in% c("gaussian", "Gamma", "inverse.gaussian",
                          "quasibinomial", "quasipoisson"))
        1 - aic/2
    else -aic/2
###    m <- if (any(x@n > 1)) x@n else w
###    sum(dbinom(round(m*origy), round(m), mu, log=TRUE))
}

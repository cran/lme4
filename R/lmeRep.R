setReplaceMethod("LMEoptimize", signature(x="lmeRep", value="list"),
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
                         ccoef(x) <- pars
                         grad <- lme4:::gradient(x, REML = value$REML, unconst = TRUE)
                         grad[constr] <- -grad[constr]/pars[constr]
                         grad
                     } else NULL
                 optimRes <- optim(st, fn, gr,
                                   method = "L-BFGS-B",
                                   lower = ifelse(constr, 1e-10, -Inf),
                                   control = list(maxit = value$msMaxIter,
                                   trace = ifelse(value$msVerbose,6,0)))
                 if (optimRes$convergence != 0) {
                     warning(paste("optim returned message",optimRes$message,"\n"))
                 }
                 ccoef(x) = optimRes$par
                 return(x)
             })

setMethod("deviance", signature(object = "lmeRep"),
          function(object, REML = FALSE, ...) {
              .Call("lmeRep_factor", object, PACKAGE = "Matrix")
              object@deviance[ifelse(REML, 2, 1)]
          })

setMethod("ranef", signature(object = "lmeRep"),
          function(object, ...) {
              .Call("lmeRep_ranef", object, PACKAGE = "Matrix")
          })

setMethod("fixef", signature(object = "lmeRep"),
          function(object, ...) {
              val = .Call("lmeRep_fixef", object, PACKAGE = "Matrix")
              names(val) = object@cnames[[".fixed"]]
              val[-length(val)]
          })

setMethod("vcov", signature(object = "lmeRep"),
          function(object, REML = TRUE, useScale = TRUE,...) {
              ## force an "lmeRep_invert"
              sc <- .Call("lmeRep_sigma", object, REML, PACKAGE = "Matrix")
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

setMethod("VarCorr", signature(x = "lmeRep"),
          function(x, REML = TRUE, useScale = TRUE, ...) {
              val = .Call("lmeRep_variances", x, PACKAGE = "Matrix")
              for (i in seq(along = val)) {
                  dimnames(val[[i]]) = list(x@cnames[[i]], x@cnames[[i]])
                  val[[i]] = as(as(val[[i]], "pdmatrix"), "corrmatrix")
              }
              new("VarCorr",
                  scale = .Call("lmeRep_sigma", x, REML),
                  reSumry = val,
                  useScale = useScale)
          })

setMethod("gradient", signature(x = "lmeRep"),
          function(x, REML, unconst, ...)
          .Call("lmeRep_gradient", x, REML, unconst))

setMethod("summary", "lmeRep",
          function(object, REML = TRUE, useScale = TRUE, ...) {
              fcoef <- fixef(object)
              corF <- as(as(vcov(object, REML, useScale), "pdmatrix"),
                         "corrmatrix")
              DF <- getFixDF(object)
              coefs <- cbind(fcoef, corF@stdDev, DF)
              nc <- object@nc
              dimnames(coefs) <-
                  list(names(fcoef), c("Estimate", "Std. Error", "DF"))
              new("summary.ssclme",
                  coefficients = as.matrix(coefs),
                  scale = .Call("lmeRep_sigma", object, REML),
                  denomDF = as.integer(DF),
                  REML = REML,
                  ngrps = unlist(lapply(object@levels, length)),
                  nobs = nc[length(nc)],
                  corFixed = corF,
                  VarCorr = VarCorr(object, REML, useScale),
                  useScale = useScale,
                  showCorrelation = FALSE)
          })

setMethod("show", "lmeRep",
          function(object) {
              fcoef <- fixef(object)
              corF <- as(as(vcov(object, REML = TRUE, useScale = TRUE), "pdmatrix"),
                         "corrmatrix")
              DF <- getFixDF(object)
              coefs <- cbind(fcoef, corF@stdDev, DF)
              nc <- object@nc
              dimnames(coefs) <-
                  list(names(fcoef), c("Estimate", "Std. Error", "DF"))
              new("summary.ssclme",
                  coefficients = as.matrix(coefs),
                  scale = .Call("lmeRep_sigma", object, REML = TRUE),
                  denomDF = as.integer(DF),
                  REML = TRUE,
                  ngrps = unlist(lapply(object@levels, length)),
                  nobs = nc[length(nc)],
                  corFixed = corF,
                  VarCorr = VarCorr(object, REML = TRUE, useScale = TRUE),
                  useScale = TRUE,
                  showCorrelation = FALSE)
          })

## calculates degrees of freedom for fixed effects Wald tests
## This is a placeholder.  The answers are generally wrong.  It will
## be very tricky to decide what a 'right' answer should be with
## crossed random effects.

setMethod("getFixDF", signature(object="lmeRep"),
          function(object, ...)
      {
          nc <- object@nc[-seq(along = object@Omega)]
          p <- nc[1] - 1
          n <- nc[2]
          rep(n-p, p)
      })

setMethod("anova", signature(object="lmeRep"),
          function(object, ...)
      {
          foo <- object
          foo@status["factored"] <- FALSE
          .Call("lmeRep_factor", foo, PACKAGE="Matrix")
          dfr <- getFixDF(foo)
          ss <- foo@RXX[ , ".response"]^2
          ssr <- ss[[".response"]]
          ss <- ss[seq(along = dfr)]
          names(ss) <- object@cnames[[".fixed"]][seq(along = dfr)]
          # FIXME: This only gives single degree of freedom tests
          ms <- ss
          df <- rep(1, length(ms))
          f <- ms/(ssr/dfr)
          P <- pf(f, df, dfr, lower.tail = FALSE)
          table <- data.frame(df, ss, ms, dfr, f, P)
          dimnames(table) <-
              list(names(ss),
                   c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
          if (any(match(names(ss), "(Intercept)", nomatch = 0)))
              table <- table[-1,]
          attr(table, "heading") <- "Analysis of Variance Table"
          class(table) <- c("anova", "data.frame")
          table
      })

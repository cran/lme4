setMethod("VarCorr", signature(x="reStruct"),
          function(x)
      {
          nobs = dim(x@original)[1]
          nfixed = length(x@random[['*fixed*']]@columns)
          denomDF <- nobs - ifelse(x@REML, nfixed, 0)
          sigma <-
              abs(x@bbetas[x@random[['*response*']]@storedRows[[1]]])/denomDF
          new("VarCorr",
              scale=sigma,
              reSumry=lapply(rev(x@random)[-c(1, 2)],
                            function(x) summary(solve(x@precision))))
      })

setMethod("VarCorr", signature(x="lme"),
          function(x)
      {
          x <- x@reStruct
          callGeneric()
      })

setMethod("VarCorr", signature(x="glmm"),
          function(x)
      {
          useScale <- !(x@family$family %in% c("binomial", "poisson"))
          x <- x@reStruct
          x <- callGeneric()
          x@useScale <- useScale
          x
      })

setMethod("show", signature(object="VarCorr"),
          function(object)
      {
          digits = max(3, getOption("digits") - 2)
          useScale = length(object@useScale) > 0 && object@useScale[1]
          sc = ifelse(useScale, object@scale,  1.)
          reStdDev <- lapply(object@reSumry, function(x, sc) sc*x@cor@stdDev,
                             sc = sc)
          reLens = unlist(lapply(reStdDev, length))
          reMat <- array('', c(sum(reLens), 4),
                         list(rep('', sum(reLens)),
                              c("Groups", "Name", "Variance", "Std.Dev.")))
          reMat[1+cumsum(reLens)-reLens, 1] = names(reLens)
          reMat[,2] = unlist(lapply(reStdDev, names))
          reMat[,3] = format(unlist(reStdDev)^2, digits = digits)
          reMat[,4] = format(unlist(reStdDev), digits = digits)
          if (any(reLens > 1) &&
              !all(sapply(object@reSumry,
                          function(x) x@noCorrelation))) {
              maxlen = max(reLens)
              corr =
                  do.call("rbind",
                          lapply(object@reSumry,
                                 function(x, maxlen) {
                                     if (x@noCorrelation) {
                                         matrix("", dim(x@cor)[1], maxlen)
                                     } else {
                                         cc = format(round(x@cor, 3),
                                                     nsmall = 3)
                                         cc[!lower.tri(cc)] = ""
                                         nr = dim(cc)[1]
                                         cbind(cc, matrix("",
                                                          nr, maxlen-nr))
                                     }
                                 }, maxlen))
              colnames(corr) = c("Corr", rep("", maxlen - 1))
              reMat = cbind(reMat, corr)
          }
          if (useScale) {
              reMat <- rbind(reMat, c("Residual", "",
                                      format(sc^2, digits = digits),
                                      format(sc, digits = digits),
                                      rep('', ncol(reMat) - 4)))
          }
          print(reMat, quote = FALSE)
      })

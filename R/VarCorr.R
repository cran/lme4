setMethod("show", signature(object="VarCorr"),
          function(object)
      {
          digits <- max(3, getOption("digits") - 2)
          useScale <- length(object@useScale) > 0 && object@useScale[1]
          sc <- ifelse(useScale, object@scale,  1.)
          reStdDev <- c(lapply(object@reSumry,
                               function(x, sc)
                               sc*x@stdDev,
                               sc = sc), list(Residual = sc))
          reLens <- unlist(c(lapply(reStdDev, length)))
          reMat <- array('', c(sum(reLens), 4),
          list(rep('', sum(reLens)),
               c("Groups", "Name", "Variance", "Std.Dev.")))
          reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
          reMat[,2] <- c(unlist(lapply(reStdDev, names)), "")
          reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
          reMat[,4] <- format(unlist(reStdDev), digits = digits)
          if (any(reLens > 1)) {
              maxlen <- max(reLens)
              corr <-
                  do.call("rbind",
                          lapply(object@reSumry,
                                 function(x, maxlen) {
                                     cc <- format(round(x, 3), nsmall = 3)
                                     cc[!lower.tri(cc)] <- ""
                                     nr <- dim(cc)[1]
                                     if (nr >= maxlen) return(cc)
                                     cbind(cc, matrix("", nr, maxlen-nr))
                                 }, maxlen))
              colnames(corr) <- c("Corr", rep("", maxlen - 1))
              reMat <- cbind(reMat, rbind(corr, rep("", ncol(corr))))
          }
          if (!useScale) reMat <- reMat[-nrow(reMat),]
          print(reMat, quote = FALSE)
      })

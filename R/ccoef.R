## constrained coefficients

## transform constrained coefficients to the Omega matrices
cco2Omega <- function(cco, k) {
    if (k < 2) return(matrix(1/cco, ncol = 1, nrow = 1))
    Lt <- diag(k)
    Lt[upper.tri(Lt)] <- cco[-(1:k)]
    ans <- t(Lt) %*% (1/(cco[1:k]) * Lt)
    ans[lower.tri(ans)] <- 0
    ans
}

## extractor
setGeneric("ccoef", function(object, ...) standardGeneric("ccoef"))

setMethod("ccoef", signature(object = "lmeRep"),
          function(object, ...) {
              nc <- object@nc
              if (length(nc) < 3) stop("inconsistent length of nc slot")
              nc <- nc[1:(length(nc)-2)]
              uco <- split(coef(object, unconst = TRUE),
                           rep(seq(a=nc), sapply(nc, function(k) (k*(k+1)/2))))
              names(uco) <- names(object@Omega)
              for (i in seq(a = nc)) {
                  uco[[i]][1:nc[i]] <- exp(-uco[[i]][1:nc[i]])
              }
              unlist(uco)
          })

setMethod("ccoef", signature(object = "ssclme"),
          function(object, ...) {
              nc <- object@nc
              if (length(nc) < 3) stop("inconsistent length of nc slot")
              nc <- nc[1:(length(nc)-2)]
              uco <- split(coef(object, unconst = TRUE),
                           rep(seq(a=nc), sapply(nc, function(k) (k*(k+1)/2))))
              names(uco) <- names(object@Omega)
              for (i in seq(a = nc)) {
                  uco[[i]][1:nc[i]] <- exp(-uco[[i]][1:nc[i]])
              }
              unlist(uco)
          })

setGeneric("ccoef<-", function(object, ..., value) standardGeneric("ccoef<-"))

setReplaceMethod("ccoef", signature(object = "lmeRep", value = "numeric"),
                 function(object, ..., value) {
                     nc <- object@nc
                     if (length(nc) < 3) stop("inconsistent length of nc slot")
                     nc <- nc[1:(length(nc)-2)]
                     Omega <- cco <-
                         split(value, rep(seq(a=nc),
                                          sapply(nc, function(k) (k*(k+1)/2))))
                     for (i in seq(a = nc)) Omega[[i]] <- cco2Omega(cco[[i]],nc[i])
                     names(Omega) <- names(object@Omega)
                     object@Omega <- Omega
                     object@status[] <- FALSE
                     object
                 })

setReplaceMethod("ccoef", signature(object = "ssclme", value = "numeric"),
                 function(object, ..., value) {
                     nc <- object@nc
                     if (length(nc) < 3) stop("inconsistent length of nc slot")
                     nc <- nc[1:(length(nc)-2)]
                     Omega <- cco <-
                         split(value, rep(seq(a=nc),
                                          sapply(nc, function(k) (k*(k+1)/2))))
                     for (i in seq(a = nc)) Omega[[i]] <- cco2Omega(cco[[i]],nc[i])
                     names(Omega) <- names(object@Omega)
                     object@Omega <- Omega
                     object@status[] <- FALSE
                     object
                 })

if (!isGeneric("lmeLevel")) {
    setGeneric("lmeLevel",
               function(precision, groups, columns, modelMatrix)
               standardGeneric("lmeLevel"))
}

setMethod("lmeLevel", signature(precision="pdMat",
                               groups="factor",
                               columns="integer",
                               modelMatrix="matrix"),
          function(precision, groups, columns, modelMatrix)
      {
          groups <- factor(groups, levels=unique(groups))
          nlev <- length(levels(groups))
          n <- length(columns)
          mat <- matrix(0.0, nrow=n, ncol=n)
          diag(mat) <- (9/64) * diag(crossprod(modelMatrix[, columns]))/nlev
          col.names <- dimnames(modelMatrix)[[2]][columns]
          dimnames(mat) <- list(col.names, col.names)
          as(precision, "matrix") <- mat
          ans <- new("lmeLevel",
                     precision=precision,
                     groups=groups,
                     columns=columns,
                     originalRows=split(as.integer(seq(length=length(groups))),
                     groups),
                     decomposedRows=vector("list", nlev),
                     storedRows=vector("list", nlev),
                     nrow=as.integer(n),
                     updateFactor=matrix(0.0, nrow=n, ncol=n),
                     nlev=as.integer(nlev))
          ans
      })

setMethod("lmeLevel", signature(precision="missing",
                               groups="missing",
                               columns="integer",
                               modelMatrix="matrix"),
          function(precision, groups, columns, modelMatrix)
      {
          ans <- new("lmeLevel",
                     precision=new("pdLogChol"),
                     groups=factor("a")[-1,drop=TRUE],
                     columns=columns,
                     originalRows=list(as.integer(seq(1, nrow(modelMatrix)))),
                     decomposedRows=vector("list", 1),
                     storedRows=vector("list", 1),
                     nrow=as.integer(0),
                     nlev=as.integer(1))
          ans
      })

setMethod("coef", signature(object="lmeLevel"),
          function(object, ...) coef(object@precision))

setReplaceMethod("coef", signature(object="lmeLevel", value="numeric"),
                 function(object, value)
             {
                 coef(object@precision) <- value
                 object
             })

setMethod("LMEhessian", signature(x="lmeLevel", A="missing",
                                  H="missing",
                                  nlev="missing"),
          function(x, A, H, nlev)
      {
          LMEhessian(x@precision, x@updateFactor, x@hessianArray,
                     nlev=x@nlev)
      })

setMethod("LMEgradient", signature(x="lmeLevel", A="missing", nlev="missing"),
          function(x, A, nlev)
      {
          LMEgradient(x@precision, x@updateFactor, x@nlev)
      })

setMethod("LMEgradient", signature(x="lmeLevel", A="matrix", nlev="missing"),
          function(x, A, nlev)
      {
          LMEgradient(x@precision, A, x@nlev)
      })

setReplaceMethod("EMupdate", signature(x="lmeLevel", nlev="missing", value="matrix"),
                 function(x, nlev, value)
             {
                 EMupdate(x@precision, x@nlev) <- value
                 x
             })

.First.lib <- function(lib, pkg) {
    if ("package:nlme" %in% search()) {
        stop(paste("Package lme4 conflicts with package nlme.\n",
                   "To attach lme4 you must restart R without package nlme."))
    }
    library.dynam(pkg, pkg, lib)
    cat(paste("This package is in development.  For production work use\n",
              "lme from package nlme or glmmPQL from package MASS.\n\n",
              "Calculation of degrees of freedom for some three-level models\n",
              "can be incorrect, as shown in example(GLMM) and in the tests\n",
              "guImmun and guPrenat.\n"))
}

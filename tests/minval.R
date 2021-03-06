if (lme4:::testLevel() > 1 || .Platform$OS.type!="windows") {
    ## example posted by Stéphane Laurent
    ## exercises bug where Nelder-Mead min objective function value was >0
    set.seed(666)
    sims <- function(I, J, sigmab0, sigmaw0){
        Mu <- rnorm(I, mean=0, sd=sigmab0)
        y <- c(sapply(Mu, function(mu) rnorm(J, mu, sigmaw0)))
        data.frame(y=y, group=gl(I,J))
    }

    I <- 3  # number of groups
    J <- 8  # number of repeats per group
    sigmab0 <- 0.15  # between standard deviation
    sigmaw0 <- 0.15  # within standard deviation

    dat <- sims(I, J, sigmab0, sigmaw0)

    library(lme4)
    isOldTol <- environment(nloptwrap)$defaultControl$xtol_abs == 1e-6

    fm3 <- lmer(y ~ (1|group), data=dat)
    stopifnot(all.equal(unname(unlist(VarCorr(fm3))),
                        switch(fm3@optinfo$optimizer,
                               "Nelder_Mead" = 0.029662844,
                               "bobyqa"      = 0.029662698,
                               "nloptwrap"   =
                                   if (isOldTol) 0.029679755 else 0.029662699,
                               stop("need new case here: value is ",unname(unlist(VarCorr(fm3))))
                               ),
                        tolerance = 1e-7))
} ## skip on windows (for speed)

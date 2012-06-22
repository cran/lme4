library(lme4)

## test AGQ=10 against glmmML results
g1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                         family = binomial, data = cbpp, nAGQ=10)
getinfo.merMod <- function(m) c(fixef(m),deviance=deviance(m))

i1 <- getinfo.merMod(g1)

if (FALSE) {
    ## comment to avoid R CMD check complaints
    ## library(glmmML)
    getinfo.glmmML <- function(m) c(coef(m),deviance=deviance(m))
    i2 <- getinfo.glmmML(glmmML(cbind(incidence, size - incidence) ~ period,
                                family = binomial,
                                cluster=herd,
                                method="ghq",
                                n.points=10,
                                data = cbpp))
}
i2 <- structure(c(-1.39924138019006, -0.991381753243723, -1.12782730029348, 
                  -1.57948092000465, 100.010029977086), .Names = c("(Intercept)", 
                                                        "period2", "period3", "period4", "deviance"))

stopifnot(all.equal(i1,i2,tol=1e-6))



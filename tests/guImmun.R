library(lme4)
data(guImmun)
numimmun <- as.numeric(guImmun$immun) - 1 # get rid of this later
system.time(fm <-
            GLMM(numimmun ~ kid2p + mom25p + ord + ethn +
                 momEd + husEd + momWork + rural + pcInd81,
                 data = guImmun, family = binomial,
                 random = list(mom = ~ 1,comm = ~1)))
summary(fm)
q("no")

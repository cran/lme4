library(lme4)
fm <- GLMM(immun ~ kid2p + mom25p + ord + ethn +
           momEd + husEd + momWork + rural + pcInd81,
           data = guImmun, family = binomial,
           random = list(mom = ~ 1,comm = ~1))
summary(fm)
q("no")

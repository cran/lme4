library(lme4)
data(guImmun)
system.time(fm <- GLMM(immun ~ kid2p + mom25p + ord + ethn +
                       momEd + husEd + momWork + rural + pcInd81,
                       data = guImmun, family = binomial,
                       random = ~ 1 |comm/mom))
summary(fm)
q("no")

library(lme4)
data(guPrenat)
system.time(fm <- GLMM(ifelse(prenat == 'Modern', 1, 0) ~ childAge +
                       motherAge + indig + momEd +
                       husEd + husEmpl + toilet + pcInd81 + ssDist,
                       data = guPrenat, family = binomial,
                       random = ~ 1|cluster/mom,
                       verbose = TRUE))
summary(fm)
q("no")


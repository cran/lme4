library(lme4)
data(s3bbx)
data(s3bby)
nres <- 1                               # change to 100 for full test
results <- array(double(nres * 11), c(nres, 11),
                 dimnames = list(NULL, c("Intercept", "chldcov", "famcov",
                 "commcov", "logsdfam", "logsdcomm", "sdfam", "sdcomm",
                 "user", "system", "elapsed")))
results[] <- NA
for (i in 1:nres) {
    try({
        thisdat <- cbind(s3bbx, modern = s3bby[,i])
        results[i, 9:11] <-
            system.time(fm <-
                        GLMM(modern ~ chldcov + famcov + commcov,
                                   data = thisdat,
                                   random = ~1|community/family,
                                   family = binomial()))[1:3]
        results[i, 1:4] <- fixef(fm)
        results[i, 5:6] <- coef(fm)
        results[i, 7:8] <- exp(results[i, 5:6])
    })
}
summary(results)
q("no")

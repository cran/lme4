library(lme4)
data(ScotsSec, package = "Matrix")
fm1 <- lme(attain ~ verbal * sex, ScotsSec,
           random = list(primary = ~1, second = ~1))
summary(fm1)
VarCorr(fm1)
fixef(fm1)
ranef(fm1)
q("no")

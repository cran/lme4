library(lme4)
fm1 <- lme(attain ~ verbal * sex, ScotsSec, ~1|primary+second,
           control = list(EMverbose=TRUE, niterEM=40, msVerbose=TRUE))
summary(fm1)
VarCorr(fm1)
fixef(fm1)
ranef(fm1)
fm2 <- lme(attain ~ verbal * sex, ScotsSec, ~1|primary+second,
           control = list(EMverbose=TRUE, msVerbose=TRUE))
summary(fm2)
q("no")

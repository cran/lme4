library(lme4)
fm1 <- lme(attain ~ verbal * sex, ScotsSec, ~1|primary+second,
           control = list(EMverbose=TRUE, niterEM=40, msVerbose=1))
summary(fm1)
VarCorr(fm1)
fixef(fm1)
ranef(fm1)
fm2 <- lme(attain ~ verbal * sex, ScotsSec, ~1|primary+second,
           control = list(EMverbose=TRUE, msVerbose=1))
summary(fm2)
lmer(attain ~ verbal * sex + (1|primary) + (1|second), ScotsSec,
     control = list(EMv = 1, msV = 1))
q("no")

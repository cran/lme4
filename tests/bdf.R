library(lme4)
fm5 <- lme(langPOST ~ IQ.ver.cen + avg.IQ.ver.cen,
                       data = bdf, random = ~ IQ.ver.cen | schoolNR)
summary(fm5)
VarCorr(fm5)
fixef(fm5)
ranef(fm5)
fm1 <- lme1(langPOST ~ IQ.ver.cen + avg.IQ.ver.cen,
                        data = bdf, random = ~ IQ.ver.cen | schoolNR)
summary(fm1)
q("no")


library(lme4)
fm5 <- lme(langPOST ~ IQ.ver.cen + avg.IQ.ver.cen,
                       data = bdf, random = ~ IQ.ver.cen | schoolNR)
summary(fm5)
VarCorr(fm5)
fixef(fm5)
ranef(fm5)
fm1 <-
    lmer(langPOST ~ IQ.ver.cen+avg.IQ.ver.cen+(IQ.ver.cen|schoolNR),
         bdf)
summary(fm1)
q("no")


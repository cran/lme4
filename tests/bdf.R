library(lme4)
data(bdf)
system.time(fm5 <- lme(langPOST ~ IQ.ver.cen + avg.IQ.ver.cen,
           data = bdf, random = ~ IQ.ver.cen | schoolNR))
summary(fm5)
fixef(fm5)
ranef(fm5)
q("no")


library(lme4)
options(show.signif.stars = FALSE)
fm1 <- lme(math ~ year, egsingle, ~1 | childid+schoolid,
           control = list(EMverbose = TRUE, msVerbose = TRUE))
summary(fm1)
fm2 <- lme(math ~ year, egsingle, list(childid=~1,schoolid=~year),
           control = list(EMverbose = TRUE, msVerbose = TRUE))
summary(fm2)
fm3 <- lme(math ~ year, egsingle, list(childid=~year,schoolid=~1),
           control = list(EMverbose = TRUE, msVerbose = TRUE))
summary(fm3)
fm4 <- lme(math ~ year, egsingle, ~year | childid+schoolid,
           control = list(EMverbose = TRUE, msVerbose = TRUE))
summary(fm4)
q("no")

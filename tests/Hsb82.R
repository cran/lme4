library(lme4)
data(Hsb82, package = "lme4")
fm1 = lme(mAch ~ meanses*cses + sector*cses,
          data = Hsb82, random = ~ cses|school)
summary(fm1)
q("no")

library(lme4)
#data(Hsb82, package = "lme4")  #Not needed in R-2.0.0 and later
lme(mAch ~ meanses*cses + sector*cses,
    data = Hsb82, random = ~ cses|school)
q("no")

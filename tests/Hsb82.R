library(lme4)
lme(mAch ~ meanses*cses + sector*cses,
    data = Hsb82, random = ~ cses|school)
lmer(mAch ~ meanses*cses + sector*cses + (cses|school), Hsb82)
q("no")

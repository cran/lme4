library(lme4)
lme(pixel ~ day + I(day^2), Pixel, list(Dog = ~ day, DS = ~ 1))
lmer(pixel ~ day + I(day^2) + (1|DS) + (day|Dog), Pixel)
q("no")

library(lme4)
Early$tos <- Early$age - 0.5
lme(cog ~ tos * trt, Early, ~ tos | id)
lmer(cog ~ tos * trt + (tos|id), Early)
q("no")

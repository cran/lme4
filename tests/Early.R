library(lme4)
Early$tos <- Early$age - 0.5
lme(cog ~ tos * trt, Early, ~ tos | id)
show(lme1(cog ~ tos * trt, Early, ~ tos | id))
q("no")

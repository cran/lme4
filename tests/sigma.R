library(lme4)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
stopifnot(all.equal(attr(VarCorr(fm1),"sc"),
		    s1 <- sigma(fm1)),
	  all.equal(s1, 25.5918, tol = 1e-6)# lme4: 25.59181783002
					#     lme4  : 25.59179822887
	  , TRUE)


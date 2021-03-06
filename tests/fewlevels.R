#### example originally from Gabor Grothendieck

source(system.file("testdata/lme-tst-funs.R", package="lme4", mustWork=TRUE))
##--> rSim.11()
testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1

if (testLevel>1) {
set.seed(1)
d1 <- rSim.11(10000, k=4)
library(nlme)
m.lme <- lme(y ~ x, random=~ 1|fac , data=d1)
(VC.lme <- VarCorr(m.lme))
detach("package:nlme")
##
library(lme4)
fm.NM <- lmer(y ~ x + (1|fac), data=d1, control=lmerControl("Nelder_Mead"))
fm.Bq <- update(fm.NM, control=lmerControl("bobyqa"))
v.lmer <- VarCorr(fm.NM)[[1]][1,1]
stopifnot(all.equal(v.lmer,19.55,tolerance=1e-3))
## was 19.5482 with old starting values (1), 19.5493 with new start algorithm
} ## testLevel>1

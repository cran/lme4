library(lme4)
set.seed(500)
Strat.mean.var.simple <- function (J, K){
   time <- rep(seq(0,1, ,length=K), J) # K measurements
   person <- rep(1:(J/2), each = K)
   treatment <- rep(0:1, each = J/2)
   treatment.ind <- rep(0:1, each = (J/2)*K)
   person1 <- rep(1:J, ,each = K)
   placebo.ind.1 <- treatment.ind < 1
   placebo.ind = ifelse( placebo.ind.1, 1, 0)
#
   mu.a.true.P = 4.8
   mu.a.true.T = 8.8
   sigma.a.true.P = 2.2
   sigma.a.true.T = 4.2
   sigma.y.true = 1.2
#
   a.true.P = rnorm (J/2, mu.a.true.P, sigma.a.true.P)
   a.true.T = rnorm (J/2, mu.a.true.T, sigma.a.true.T)
#
   y.P <- rnorm( (J/2)*K, a.true.P[person], sigma.y.true)
   y.T <- rnorm( (J/2)*K, a.true.T[person], sigma.y.true)
   y <- c(y.P, y.T)
   return ( data.frame (y, time, person1, treatment.ind, placebo.ind))
}
fake <- Strat.mean.var.simple (J = 10, K = 4)
lme.no.strat <- lmer(y ~  treatment.ind + (1 | person1) , data = fake)
lme.strat <- lmer (y ~  treatment.ind + ( 0 + placebo.ind | person1) + (0 +
treatment.ind | person1), data = fake)
try(coef(lme.strat))
coef(lme.no.strat)



## ----preliminaries,include=FALSE,cache=FALSE,message=FALSE----------
options(width=70, show.signif.stars=FALSE,
        str=strOptions(strict.width="cut"),
        ## prefer empty continuation for reader's cut'n'paste:
        continue = "   ", #JSS: prompt = "R> ", continue = "+  ",
        useFancyQuotes = FALSE)
library("knitr")
library("lme4")
library("ggplot2")# Keeping default theme, nicer "on line":
#JSS theme_set(theme_bw())
library("grid")
zmargin <- theme(panel.spacing=unit(0,"lines"))
library("lattice")
library("minqa")
library("reformulas")
opts_chunk$set(engine='R',dev='pdf', fig.width=9, fig.height=5.5,
               prompt=TRUE, cache=TRUE, tidy=FALSE, comment=NA, error = FALSE)
knitr::render_sweave()

## ----phiToTheta, cache = FALSE--------------------------------------
phiToTheta <- function(phi) {
    theta5 <- -(phi[2]*phi[3])/phi[4]
    c(phi[1:4], theta5, phi[5])
}

## ----compute deviance function modular, cache = FALSE---------------
lf <- lFormula(formula = form, data = dat, REML = TRUE)
devf <- do.call(mkLmerDevfun, lf)

## ----wrapper modular, cache = FALSE---------------------------------
devfWrap <- function(phi) devf(phiToTheta(phi))

## ----opt modular, cache = FALSE-------------------------------------
opt <- bobyqa(par = lf$reTrms$theta[-5],
  fn = devfWrap,
  lower = lf$reTrms$lower[-5])

## ----varCorr modular, cache = FALSE---------------------------------
vcEst <- vec2mlist(Cv_to_Vv(phiToTheta(opt$par)))[[1]]
dimnames(vcEst) <- rep(lf$reTrms$cnms, 2)
round(vcEst, 2)
vc


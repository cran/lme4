#ifndef LME4_LMER2_H
#define LME4_LMER2_H

#include "lme4_utils.h"
#include <R_ext/stats_package.h>

SEXP glmer_bhat(SEXP x);
SEXP glmer_bhat2(SEXP x);
SEXP glmer_eta(SEXP x);
SEXP glmer_reweight(SEXP x);

SEXP lmer2_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp,
		   SEXP verbose, SEXP deviance);
SEXP lmer2_create(SEXP fr, SEXP FL, SEXP Ztl, SEXP glmp,
		  SEXP method, SEXP mc, SEXP mod);
SEXP lmer2_deviance(SEXP x, SEXP which);
SEXP lmer2_getPars(SEXP x);
SEXP lmer2_optimize(SEXP x, SEXP verb);
SEXP lmer2_postVar(SEXP x);
SEXP lmer2_ranef(SEXP x);
SEXP lmer2_setPars(SEXP x, SEXP pars);
SEXP lmer2_sigma(SEXP x, SEXP which);
SEXP lmer2_update_effects(SEXP x);
SEXP lmer2_update_y(SEXP x, SEXP yp);
SEXP lmer2_validate(SEXP x);
SEXP lmer2_vcov(SEXP x);

SEXP nlmer_bhat(SEXP x);
SEXP nlmer_create(SEXP env, SEXP model, SEXP frame, SEXP pnames,
		  SEXP call, SEXP flist, SEXP Xt, SEXP Zt, SEXP y,
		  SEXP weights, SEXP cnames, SEXP Gp, SEXP fixef);
SEXP nlmer_eval_model(SEXP x, SEXP uform);
SEXP nlmer_optimize(SEXP x, SEXP verb);
SEXP nlmer_update_Mt(SEXP x);
SEXP nlmer_update_Vt(SEXP x);
SEXP nlmer_update_ranef(SEXP x);
SEXP nlmer_update_wrkres(SEXP x);
SEXP nlmer_validate(SEXP x);

#endif /* LME4_LMER2_H */

#ifndef LME4_LMER_H
#define LME4_LMER_H

#include "lme4_utils.h"
#include <Rmath.h>
#include <Rversion.h>
#include <R_ext/Lapack.h>
#include <R_ext/Constants.h>
#include <R_ext/Random.h>

SEXP glmer_MCMCsamp(SEXP GSp, SEXP savebp, SEXP nsampp, SEXP transp, SEXP verbose);
SEXP glmer_PQL(SEXP GSp);
SEXP glmer_devLaplace(SEXP pars, SEXP GSp);
SEXP glmer_finalize(SEXP GSpt);
SEXP glmer_init(SEXP rho);

SEXP lme4_rWishart(SEXP ns, SEXP df, SEXP scal);

SEXP mer_ECMEsteps(SEXP x, SEXP nsteps, SEXP Verbp);
/* SEXP mer_Hessian(SEXP x);  not yet */
SEXP mer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp, SEXP verbose);
SEXP mer_coef(SEXP x, SEXP pType);
SEXP mer_coefGets(SEXP x, SEXP coef, SEXP pType);
SEXP mer_create(SEXP fl, SEXP Zt, SEXP Xp, SEXP yp, SEXP REMLp, SEXP nc, SEXP cnames);
SEXP mer_dtCMatrix(SEXP x);
SEXP mer_dtCMatrix_inv(SEXP x);
SEXP mer_factor(SEXP x);
SEXP mer_fitted(SEXP x);
SEXP mer_fixef(SEXP x);
SEXP mer_gradComp(SEXP x);
SEXP mer_gradient(SEXP x, SEXP pType);
SEXP mer_hat_trace(SEXP x);
SEXP mer_hat_trace2(SEXP x);
SEXP mer_initial(SEXP x);
SEXP mer_isNested(SEXP x);
SEXP mer_postVar(SEXP x);
SEXP mer_ranef(SEXP x);
SEXP mer_secondary(SEXP x);
SEXP mer_sigma(SEXP x, SEXP REML);
SEXP mer_simulate(SEXP x, SEXP nsimP);
SEXP mer_update_ZXy(SEXP x);
SEXP mer_update_y(SEXP x, SEXP ynew);
SEXP mer_validate(SEXP x);

/* SEXP Zt_create(SEXP fl, SEXP Ztl); */
/* SEXP Zt_create1(SEXP fl, SEXP Ztl); */
SEXP Ztl_sparse(SEXP fl, SEXP Ztl);
SEXP Zt_carryOver(SEXP f, SEXP Zt);

#endif /* LME4_LMER_H */

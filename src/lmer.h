#ifndef LME4_LMER_H
#define LME4_LMER_H

#include "lme4_utils.h"
#include "Wishart.h"

SEXP mer_ECMEsteps(SEXP x, SEXP nsteps, SEXP Verbp);
/* SEXP mer_Hessian(SEXP x);  not yet */
SEXP mer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp, SEXP verbose);
SEXP mer_coef(SEXP x, SEXP pType);
SEXP mer_coefGets(SEXP x, SEXP coef, SEXP pType);
SEXP mer_create(SEXP fl, SEXP Zt, SEXP Xp, SEXP yp, SEXP REMLp, SEXP nc, SEXP cnames);
SEXP mer_dtCMatrix(SEXP x);
SEXP mer_dtCMatrix_inv(SEXP x);
SEXP mer_fitted(SEXP x);
SEXP mer_fixef(SEXP x);
SEXP mer_gradient(SEXP x, SEXP pType);
SEXP mer_hat_trace(SEXP x);
SEXP mer_hat_trace2(SEXP x);
SEXP mer_initial(SEXP x);
SEXP mer_isNested(SEXP x);
SEXP mer_postVar(SEXP x);
SEXP mer_ranef(SEXP x);
SEXP mer_sigma(SEXP x, SEXP REML);
SEXP mer_simulate(SEXP x, SEXP nsimP);
SEXP mer_update_ZXy(SEXP x);
SEXP mer_update_y(SEXP x, SEXP ynew);
SEXP mer_validate(SEXP x);

SEXP Ztl_sparse(SEXP fl, SEXP Ztl);
SEXP Zt_carryOver(SEXP f, SEXP Zt, SEXP tvar, SEXP discount);

#endif /* LME4_LMER_H */

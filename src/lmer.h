#ifndef LME4_LMER_H
#define LME4_LMER_H

#include "lme4_utils.h"
#include "Wishart.h"
				/* positions in the deviance vector */
extern const char attr_hidden *DEVIANCE_NAMES[];
			 /* {"ML", "REML", "ldZ", "ldX", "lr2", ""} */
#define ML_POS   0
#define REML_POS 1
#define ldZ_POS  2
#define ldX_POS  3
#define lr2_POS  4
				/* positions in the dims vector */
extern const char attr_hidden *DIMS_NAMES[];
	      /* {"nf", "n", "p", "q", "REML", "glmm", ""} */
#define nf_POS     0
#define n_POS      1
#define p_POS      2
#define q_POS      3
/* FIXME: change the next two codes to a single code "type" with
   values type < 0 => LMM & REML, type == 0 => LMM & ML,
   type == 1 => GLMM & Laplace, type == 2 => GLMM & AGQ
*/
#define isREML_POS 4
#define glmm_POS 5

#define isREML(x) INTEGER(GET_SLOT(x, lme4_dimsSym))[isREML_POS]
#define isGLMM(x) INTEGER(GET_SLOT(x, lme4_dimsSym))[glmm_POS]

SEXP mer_ECMEsteps(SEXP x, SEXP nsteps, SEXP Verbp);
SEXP mer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp,
		  SEXP verbose, SEXP deviance);
SEXP mer_coef(SEXP x, SEXP pType);
SEXP mer_coefGets(SEXP x, SEXP coef, SEXP pType);
SEXP mer_create(SEXP fl, SEXP Zt, SEXP Xp, SEXP yp, SEXP REMLp, SEXP nc, SEXP cnames);
SEXP mer2_create(SEXP fl, SEXP Zt, SEXP Xp, SEXP yp, SEXP REMLp,
		SEXP nc, SEXP cnames, SEXP offset, SEXP wts);
SEXP mer_dtCMatrix(SEXP x);
SEXP mer_dtCMatrix_inv(SEXP x);
SEXP mer2_deviance(SEXP x, SEXP which);
SEXP mer_fitted(SEXP x);
SEXP mer_fixef(SEXP x);
SEXP mer2_getPars(SEXP x);
SEXP mer_gradient(SEXP x, SEXP pType);
SEXP mer_hat_trace(SEXP x);
SEXP mer_hat_trace2(SEXP x);
SEXP mer_initial(SEXP x);
SEXP mer_isNested(SEXP x);
SEXP mer_postVar(SEXP x);
SEXP mer_ranef(SEXP x);
SEXP mer2_setPars(SEXP x, SEXP pars);
SEXP mer_sigma(SEXP x, SEXP REML);
SEXP mer2_sigma(SEXP x, SEXP which);
SEXP mer_simulate(SEXP x, SEXP nsimP);
SEXP mer2_update_effects(SEXP x);
SEXP mer_update_ZXy(SEXP x);
SEXP mer_update_y(SEXP x, SEXP ynew);
SEXP mer_validate(SEXP x);
SEXP mer2_vcov(SEXP x);

SEXP Ztl_sparse(SEXP fl, SEXP Ztl);
SEXP Zt_carryOver(SEXP f, SEXP Zt, SEXP tvar, SEXP discount);

#endif /* LME4_LMER_H */

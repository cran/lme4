#include "lme4_utils.h"
#include "glmer.h"
#include "lmer2.h"
#include "pedigree.h"
#include <R_ext/Rdynload.h>

/* Syms.h needs to be included a second time (it is already included
 * through lme4_utils.h) so the symbols are defined without extern.
 */
#include "Syms.h" 

static R_CallMethodDef CallEntries[] = {
    {"glmer_MCMCsamp", (DL_FUNC) &glmer_MCMCsamp, 6},
    {"glmer_PQL", (DL_FUNC) &glmer_PQL, 1},
    {"glmer_bhat", (DL_FUNC) &glmer_bhat, 1},
    {"glmer_bhat2", (DL_FUNC) &glmer_bhat2, 1},
    {"glmer_devLaplace", (DL_FUNC) &glmer_devLaplace, 2},
    {"glmer_finalize", (DL_FUNC) &glmer_finalize, 1},
    {"glmer_init", (DL_FUNC) &glmer_init, 2},
    {"glmer_eta", (DL_FUNC) &glmer_eta, 1},
    {"glmer_reweight", (DL_FUNC) &glmer_reweight, 1},

    {"lme4_rWishart", (DL_FUNC) &lme4_rWishart, 3},
    {"lmer2_create", (DL_FUNC) &lmer2_create, 7},
    {"lmer2_MCMCsamp", (DL_FUNC) &lmer2_MCMCsamp, 6},
    {"lmer2_deviance", (DL_FUNC) &lmer2_deviance, 2},
    {"lmer2_getPars", (DL_FUNC) &lmer2_getPars, 1},
    {"lmer2_optimize", (DL_FUNC) &lmer2_optimize, 2},
    {"lmer2_postVar", (DL_FUNC) &lmer2_postVar, 1},
    {"lmer2_ranef", (DL_FUNC) &lmer2_ranef, 1},
    {"lmer2_setPars", (DL_FUNC) &lmer2_setPars, 2},
    {"lmer2_sigma", (DL_FUNC) &lmer2_sigma, 2},
    {"lmer2_update_effects", (DL_FUNC) &lmer2_update_effects, 1},
    {"lmer2_update_y", (DL_FUNC) &lmer2_update_y, 2},
    {"lmer2_validate", (DL_FUNC) &lmer2_validate, 1},
    {"lmer2_vcov", (DL_FUNC) &lmer2_vcov, 1},

    {"mer_ECMEsteps", (DL_FUNC) &mer_ECMEsteps, 3},
    {"mer_MCMCsamp", (DL_FUNC) &mer_MCMCsamp, 6},
    {"mer_coef", (DL_FUNC) &mer_coef, 2},
    {"mer_coefGets", (DL_FUNC) &mer_coefGets, 3},
    {"mer_create", (DL_FUNC) &mer_create, 7},
    {"mer_dtCMatrix", (DL_FUNC) &mer_dtCMatrix, 1},
    {"mer_dtCMatrix_inv", (DL_FUNC) &mer_dtCMatrix_inv, 1},
    {"mer_factor", (DL_FUNC) &mer_factor, 1},
    {"mer_fitted", (DL_FUNC) &mer_fitted, 1},
    {"mer_fixef", (DL_FUNC) &mer_fixef, 1},
    {"mer_gradComp", (DL_FUNC) &mer_gradComp, 1},
    {"mer_gradient", (DL_FUNC) &mer_gradient, 2},
    {"mer_hat_trace", (DL_FUNC) &mer_hat_trace, 1},
    {"mer_hat_trace2", (DL_FUNC) &mer_hat_trace2, 1},
    {"mer_initial", (DL_FUNC) &mer_initial, 1},
    {"mer_isNested", (DL_FUNC) &mer_isNested, 1},
    {"mer_postVar", (DL_FUNC) &mer_postVar, 1},
    {"mer_ranef", (DL_FUNC) &mer_ranef, 1},
    {"mer_secondary", (DL_FUNC) &mer_secondary, 1},
    {"mer_sigma", (DL_FUNC) &mer_sigma, 2},
    {"mer_simulate", (DL_FUNC) &mer_simulate, 2},
    {"mer_update_ZXy", (DL_FUNC) &mer_update_ZXy, 1},
    {"mer_update_y", (DL_FUNC) &mer_update_y, 2},
    {"mer_validate", (DL_FUNC) &mer_validate, 1},

    {"nlmer_bhat", (DL_FUNC) &nlmer_bhat, 1},
    {"nlmer_create", (DL_FUNC) &nlmer_create, 13},
    {"nlmer_eval_model", (DL_FUNC) &nlmer_eval_model, 2},
    {"nlmer_optimize", (DL_FUNC) &nlmer_optimize, 2},
    {"nlmer_update_Mt", (DL_FUNC) &nlmer_update_Mt, 1},
    {"nlmer_update_Vt", (DL_FUNC) &nlmer_update_Vt, 1},
    {"nlmer_update_ranef", (DL_FUNC) &nlmer_update_ranef, 1},
    {"nlmer_update_wrkres", (DL_FUNC) &nlmer_update_wrkres, 1},
    {"nlmer_validate", (DL_FUNC) &nlmer_validate, 1},

    {"pedigree_chol", (DL_FUNC) &pedigree_chol, 2},

    {"Zt_carryOver", (DL_FUNC) &Zt_carryOver, 4},
    {"Ztl_sparse", (DL_FUNC) &Ztl_sparse, 2},

    {NULL, NULL, 0}
};

cholmod_common c;

#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
void R_init_lme4(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);


    M_R_cholmod_start(&c);
    c.final_ll = 0;	    /* LDL form of simplicial factorization */

    lme4_ASym = install("A");
    lme4_DSym = install("D");
    lme4_DimSym = install("Dim");
    lme4_DimNamesSym = install("Dimnames");
    lme4_GpSym = install("Gp");
    lme4_LSym = install("L");
    lme4_OmegaSym = install("Omega");
    lme4_RXXSym = install("RXX");
    lme4_RZXSym = install("RZX");
    lme4_RZXinvSym = install("RZXinv");
    lme4_STSym = install("ST");
    lme4_XSym = install("X");
    lme4_XtXSym = install("XtX");
    lme4_XtySym = install("Xty");
    lme4_ZZpOSym = install("ZZpO");
    lme4_ZtSym = install("Zt");
    lme4_ZtXSym = install("ZtX");
    lme4_ZtZSym = install("ZtZ");
    lme4_ZtySym = install("Zty");
    lme4_ZXytSym = install("ZXyt");
    lme4_bVarSym = install("bVar");
    lme4_cnamesSym = install("cnames");
    lme4_devCompSym = install("devComp");
    lme4_devianceSym = install("deviance");
    lme4_diagSym = install("diag");
    lme4_dimsSym = install("dims");
    lme4_etaSym = install("eta");
    lme4_factorSym = install("factor");
    lme4_fixefSym = install("fixef");
    lme4_flistSym = install("flist");
    lme4_gradCompSym = install("gradComp");
    lme4_iSym = install("i");
    lme4_muSym = install("mu");
    lme4_ncSym = install("nc");
    lme4_offsetSym = install("offset");
    lme4_pSym = install("p");
    lme4_permSym = install("perm");
    lme4_rXySym = install("rXy");
    lme4_rZySym = install("rZy");
    lme4_ranefSym = install("ranef");
    lme4_statusSym = install("status");
    lme4_uploSym = install("uplo");
    lme4_weightsSym = install("weights");
    lme4_wrkresSym = install("wrkres");
    lme4_wtsSym = install("wts");
    lme4_xSym = install("x");
    lme4_ySym = install("y");
}

void R_unload_lme4(DllInfo *dll){
    M_cholmod_finish(&c);
}

#include "lmer.h"
#include <R_ext/Rdynload.h>
#include "Matrix.h"
#include "Syms.h" 

static R_CallMethodDef CallEntries[] = {
/*     {"lme4_rWishart", (DL_FUNC) &lme4_rWishart, 3}, */

    {"mer_MCMCsamp", (DL_FUNC) &mer_MCMCsamp, 2},
    {"mer_ST_chol", (DL_FUNC) &mer_ST_chol, 1},
    {"mer_ST_getPars", (DL_FUNC) &mer_ST_getPars, 1},
    {"mer_ST_initialize", (DL_FUNC) &mer_ST_initialize, 3},
    {"mer_ST_setPars", (DL_FUNC) &mer_ST_setPars, 2},
    {"mer_create_L", (DL_FUNC) &mer_create_L, 1},
    {"mer_optimize", (DL_FUNC) &mer_optimize, 2},
    {"mer_postVar", (DL_FUNC) &mer_postVar, 1},
    {"mer_update_L", (DL_FUNC) &mer_update_L, 1},
    {"mer_update_RX", (DL_FUNC) &mer_update_RX, 1},
    {"mer_update_dev", (DL_FUNC) &mer_update_dev, 1},
    {"mer_update_ranef", (DL_FUNC) &mer_update_ranef, 1},
    {"mer_update_mu", (DL_FUNC) &mer_update_mu, 1},
    {"mer_update_u", (DL_FUNC) &mer_update_u, 2},
    {"mer_validate", (DL_FUNC) &mer_validate, 1},

    {"merMCMC_validate", (DL_FUNC) &merMCMC_validate, 1},
    {"merMCMC_VarCorr", (DL_FUNC) &merMCMC_VarCorr, 2},

    {"spR_optimize", (DL_FUNC) &spR_optimize, 2},
    {"spR_update_mu", (DL_FUNC) &spR_update_mu, 1},
/*     {"Zt_carryOver", (DL_FUNC) &Zt_carryOver, 4}, */

    {NULL, NULL, 0}
};

/** cholmod_common struct local to the package */
cholmod_common c;

/** Initializer for lme4, called upon loading the package.
 *
 *  Register routines that can be called directly from R.
 *  Initialize CHOLMOD and require the LL' form of the factorization.
 *  Install the symbols to be used by functions in the package.
 */
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
void R_init_lme4(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);


    M_R_cholmod_start(&c);
    c.final_ll = 1;	    /* LL' form of simplicial factorization */

    lme4_ASym = install("A");
    lme4_CmSym = install("Cm");
    lme4_CxSym = install("Cx");
    lme4_DimSym = install("Dim");
    lme4_GpSym = install("Gp");
    lme4_LSym = install("L");
    lme4_RXSym = install("RX");
    lme4_RZXSym = install("RZX");
    lme4_STSym = install("ST");
    lme4_VSym = install("V");
    lme4_XSym = install("X");
    lme4_ZtSym = install("Zt");
    lme4_devianceSym = install("deviance");
    lme4_dimsSym = install("dims");
    lme4_envSym = install("env");
    lme4_etaSym = install("eta");
    lme4_fixefSym = install("fixef");
    lme4_flistSym = install("flist");
    lme4_gradientSym = install("gradient");
    lme4_iSym = install("i");
    lme4_ncSym = install("nc");
    lme4_nlmodelSym = install("nlmodel");
    lme4_muEtaSym = install("muEta");
    lme4_muSym = install("mu");
    lme4_offsetSym = install("offset");
    lme4_pSym = install("p");
    lme4_permSym = install("perm");
    lme4_pWtSym = install("pWt");
    lme4_ranefSym = install("ranef");
    lme4_residSym = install("resid");
    lme4_sigmaSym = install("sigma");
    lme4_sqrtrWtSym = install("sqrtrWt");
    lme4_sqrtXWtSym = install("sqrtXWt");
    lme4_uSym = install("u");
    lme4_varSym = install("var");
    lme4_xSym = install("x");
    lme4_ySym = install("y");
}

/** Finalizer for lme4 called upon unloading the package.
 *
 */
void R_unload_lme4(DllInfo *dll){
    M_cholmod_finish(&c);
}

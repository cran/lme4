/*
   Add the foolowing files from the PORT library to enable
   support for the ms optimizer.

   da7sst.f dd7dog.f dd7tpr.f ditsum.f divset.f dl7itv.f dl7ivm.f
   dl7tvm.f dl7upd.f dl7vml.f dmnf.f dmng.f dparck.f dr7mdc.f
   drldst.f drmnf.f drmng.f ds7grd.f dv2axy.f dv2nrm.f dv7cpy.f
   dv7dfl.f dv7scp.f dv7vmp.f dw7zbf.f i7mdcn.f stopx.f

*/

#include "reStruct.h"

/* Uncomment the following to add support for the Port library 
   optimizer. The appropriate FORTRAN files from the PORT library
   must also be compiled in to do this. */

/* #define NLME_USE_PORTLIB */

#ifdef NLME_USE_PORTLIB

void F77_NAME(dmng)(int* N, double* D, double* X,
                    void (*CALCF)(int*, double*, int*, double*,
                                  int*, double*, void (*ufparm)()),
                    void (*CALCG)(int*, double*, int*, double*,
                                  int*, double*, void (*ufparm)()),
                    int* IV, int* LIV, int* LV, double* V,
                    int* UIPARM, double* URPARM, void (*UFPARM)());

void F77_NAME(dmnf)(int* N, double* D, double* X,
                    void (*CALCF)(int*, double*, int*, double*,
                                  int*, double*, void (*ufparm)()),
                    int* IV, int* LIV, int* LV, double* V,
                    int* UIPARM, double* URPARM, void (*UFPARM)());

static void
nlme_calcf(int* n, double* x, int* nf, double* f, int* uiparm,
           double* urparm, void (*UFPARM)())
{
    SEXP pars = (SEXP) uiparm;
    SEXP reStruct = (SEXP) urparm;
    if (REAL(pars) != x) {
        int i;
        for (i = 0; i < *n; i++)
            REAL(pars)[i] = x[i];
    }
    nlme_setParameters(&reStruct, pars);
    *f = -nlme_logLikelihood_internal(reStruct, 0, 0);
}

static void
nlme_calcg(int* n, double* x, int* nf, double* g, int* uiparm,
           double* urparm, void (*UFPARM)())
{
    SEXP pars = (SEXP) uiparm;
    SEXP reStruct = (SEXP) urparm;
    double* grad;
    int i;
    if (REAL(pars) != x) {
        for (i = 0; i < *n; i++)
            REAL(pars)[i] = x[i];
    }
    nlme_setParameters(&reStruct, pars);
    grad = REAL(eval(lang2(install("LMEgradient"), reStruct), R_GlobalEnv));
    for (i = 0; i < *n; i++) {
        g[i] = -grad[i];
    }
}

SEXP
nlme_msOptimize(SEXP msMaxIter, SEXP msTol, SEXP scale, SEXP msVerbose,
                SEXP reStruct, SEXP pars, SEXP useAnalyticGrad)
{
    int analyticGrad = asLogical(useAnalyticGrad);
    int two = 2;
    int npars = LENGTH(pars);
    double* d = REAL(scale);
    double* x = REAL(pars);
    int liv = 60;
    int *iv = Calloc(liv, int);
    int lv = analyticGrad?(71+npars*(npars+15)/2):(77+npars*(npars+17)/2);
    double *v = Calloc(lv, double);
    int* uiparm = (int*) pars;
    double* urparm;
    if (NAMED(reStruct) && !asLogical(GET_SLOT(reStruct,
                                               install("dontCopy")))) {
        reStruct = duplicate(reStruct);
    }
    urparm = (double*) PROTECT(reStruct);
    F77_CALL(divset)(&two, iv, &liv, &lv, v);
    iv[17] = asInteger(msMaxIter);
    iv[21] = 0;
    iv[22] = -1;
    iv[23] = 0;
    v[25] = 0.001;
    v[30] = 0.0;
    v[33] = asReal(msTol);
    if (analyticGrad)
        F77_CALL(dmng)(&npars, d, x, nlme_calcf, nlme_calcg, iv, &liv, &lv, v,
                       uiparm, urparm, NULL);
    else F77_CALL(dmnf)(&npars, d, x, nlme_calcf, iv, &liv, &lv, v,
                        uiparm, urparm, NULL);
    liv = iv[0];
    Free(iv);
    Free(v);
    if (liv > 11) {
        error("problem in ms");
    }
    switch(liv) {
    case 7:
        warning("singular convergence in ms");
        break;
    case 8:
        warning("false convergence in ms");
        break;
    case 9:
        warning("function evaluation limit reached without convergence");
        break;
    case 10:
        warning("iteration limit reached without convergence in ms");
        break;
    case 11:
        warning("external interrupt to ms");
        break;
    default:
        break;
    }
    reStruct = nlme_commonDecompose(reStruct, pars);
    UNPROTECT(1);
    return reStruct;
}

#endif

#ifndef LME4_GLMER_H
#define LME4_GLMER_H

#include "lmer.h"

typedef struct glmer_struct
{
    SEXP cv;         /* control values */
    SEXP mer;	     /* mixed-effects representation */
    SEXP rho;        /* environment in which to evaluate the calls */
    SEXP eta;        /* linear predictor */
    SEXP mu;         /* mean vector */
    SEXP LMEopt;     /* expression for LME optimization */
    SEXP dev_resfunc; /* expression for deviance residuals (if fltype == 0) */
    SEXP linkinv;    /* expression for inverse link evaluation (if fltype == 0) */
    SEXP mu_eta;     /* expression for dmu/deta evaluation (if fltype == 0) */
    SEXP vfunc;      /* expression for variance evaluation (if fltype == 0) */
    double *dev_res; /* deviance residuals */
    double *dmu_deta;/* derivative vector */
    double *var;     /* variance vector */
    double *offset;  /* offset for GLM */
    double *wts;     /* prior weights for GLM */
    double *y;       /* copy of response vector */
    double *etaold;  /* previous value of eta for evaluating  */
    int fltype;      /* family-link type */
    int n;	     /* length of the response vector */
    int p;	     /* length of fixed effects vector */
    int nf;	     /* number of grouping factors */
    int npar;        /* total length of the parameter */
    int niterEM;     /* default number of ECME iterations */
    int EMverbose;   /* logical indicator */
    int maxiter;     /* maximum number of IRLS iterations */
    double tol;      /* convergence tolerance for IRLS iterations */
} glmer_struct, *GlmerStruct;

SEXP glmer_MCMCsamp(SEXP GSp, SEXP savebp, SEXP nsampp, SEXP transp,
		    SEXP verbose, SEXP deviancep);
SEXP glmer_PQL(SEXP GSp);
SEXP glmer_devLaplace(SEXP pars, SEXP GSp);
SEXP glmer_finalize(SEXP GSpt);
SEXP glmer_init(SEXP rho, SEXP fltype);

#endif /* LME4_LMER_H */

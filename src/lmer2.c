#include "lmer2.h"

/* TODO: */
/* Consider the steps in reimplementing AGQ.  First you need to find
   bhat, then evaluate the posterior variances, then step out
   according to the posterior variance, evaluate the integrand
   relative to the step. */
/* Because the Gauss-Hermite quadrature is formed as a sum, it is
 * necessary to divide the contributions to the deviance according to
 * the levels of the random effects.  This means that it is only
 * practical to use AGQ when the response vector can be split into
 * sections that are conditionally independent. As far as I can see
 * this will mean a single grouping factor only. */

/* Functions for the lmer2 representation */

				/* positions in the deviance vector */
enum devP {ML_POS=0, REML_POS, ldZ_POS, ldX_POS, lr2_POS, bqd_POS, Sdr_POS};
			/* {"ML", "REML", "ldZ", "ldX", "lr2", "bQuad" "sumDevR" ""} */
				/* positions in the dims vector */
enum dimP {nf_POS=0, n_POS, p_POS, q_POS, s_POS, np_POS, isREML_POS, famType_POS, isNest_POS};
	      /* {"nf", "n", "p", "q", "s", "np", "isREML", "famType", "isNested"} */

#define isREML(x) INTEGER(GET_SLOT(x, lme4_dimsSym))[isREML_POS]
#define isGLMM(x) (INTEGER(GET_SLOT(x, lme4_dimsSym))[famType_POS] >= 0)
#define isNested(x) INTEGER(GET_SLOT(x, lme4_dimsSym))[isNest_POS]

/**
 * Return the element of a given name from a named list
 *
 * @param list
 * @param nm name of desired element
 *
 * @return element of list with name nm
 */
static SEXP R_INLINE getListElement(SEXP list, char *nm) {
    SEXP names = getAttrib(list, R_NamesSymbol);
    int i;

    if (!isNull(names))
	for (i = 0; i < LENGTH(names); i++)
	    if (!strcmp(CHAR(STRING_ELT(names, i)), nm))
		return(VECTOR_ELT(list, i));
    return R_NilValue;
}

/**
 * Check the ZtZ matrix to see if it is a simple design from a nested
 * sequence of grouping factors.
 *
 * @param nf number of factors
 * @param nc[] number of columns per factor
 * @param Gp[] group pointers
 * @param p[] column pointers for the lower triangle of ZtZ
 *
 * @return 1 for a simple nested sequence, 0 otherwise.
 */
static int check_nesting(int nf, SEXP ST, const int Gp[], const int p[])
{
    int **cnz = Calloc(nf, int*), *nc = Calloc(nf, int), ans = 1, i, j, k, nct;

    for (i = 0, nct = 0; i < nf; i++) { /* total number of columns */
	nc[i] = *INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol));
	nct += nc[i];
	cnz[i] = Calloc(nc[i], int);
    }
    for (i = 0; i < nf; i++) {	/* target number of nonzeros per column */
	for (j = 0; j < nc[i]; j++) cnz[i][j] = nct - j;
	nct -= nc[i];
    }
    for (i = 0; i < nf && ans; i++) { /* check for consistent nonzeros*/
	int nlev = (Gp[i + 1] - Gp[i])/nc[i];
	for (j = 0; j < nlev && ans; j++) {
	    for (k = 0; k < nc[i] && ans; k++) {
		int jj =  Gp[i] + j * nc[i] + k; /* column in ZtZ */
		if ((p[jj + 1] - p[jj]) != cnz[i][k]) ans = 0;
	    }
	}
    }
    for (i = 0, nct = 0; i < nf; i++) Free(cnz[i]);
    Free(cnz); Free(nc);
    return ans;
}

				/* Utilities for generalized linear mixed models */
static const double LTHRESH = 30.;
static const double MLTHRESH = -30.;
static double MPTHRESH = 0;
static double PTHRESH = 0;
static const double INVEPS = 1/DOUBLE_EPS;

/** 
 * Evaluate x/(1 - x). An inline function is used so that x is
 * evaluated once only. 
 * 
 * @param x input in the range (0, 1)
 * 
 * @return x/(1 - x) 
 */
static R_INLINE double x_d_omx(double x) {
    if (x < 0 || x > 1)
	error(_("Value %d out of range (0, 1)"), x);
    return x/(1 - x);
}

/** 
 * Evaluate x/(1 + x). An inline function is used so that x is
 * evaluated once only.
 * 
 * @param x input
 * 
 * @return x/(1 + x) 
 */
static R_INLINE double x_d_opx(double x) {return x/(1 + x);}

static R_INLINE double y_log_y(double y, double mu)
{
    return (y) ? (y * log(y/mu)) : 0;
}

/**
 * Evaluate the inverse link function at eta storing the result in mu
 *
 * @param x pointer to a glmer2 object
 */
static void glmer_linkinv(SEXP x)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, n = dims[n_POS], fltype = dims[famType_POS];
    double *eta = REAL(GET_SLOT(x, lme4_etaSym)),
	*mu = REAL(GET_SLOT(x, lme4_muSym));

    switch(fltype) {
    case 1: 			/* binomial with logit link */
	for (i = 0; i < n; i++) {
	    double etai = eta[i], tmp;
	    tmp = (etai < MLTHRESH) ? DOUBLE_EPS :
		((etai > LTHRESH) ? INVEPS : exp(etai));
	    mu[i] = x_d_opx(tmp);
	}
	break;
    case 2:			/* binomial with probit link */
	if (!MPTHRESH) {
	    MPTHRESH = qnorm5(DOUBLE_EPS, 0, 1, 1, 0);
	    PTHRESH = -MPTHRESH;
	}
	for (i = 0; i < n; i++) {
	    double etai = eta[i];
	    mu[i] = (etai < MPTHRESH) ? DOUBLE_EPS :
		((etai > PTHRESH) ? 1 - DOUBLE_EPS :
		 pnorm5(etai, 0, 1, 1, 0));
	}
	break;
    case 3:			/* Poisson with log link */
	for (i = 0; i < n; i++) {
	    double tmp = exp(eta[i]);
	    mu[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
	}
 	break;
    default:
	error(_("General form of glmer_linkinv not yet written"));
     } 
}

/**
 * Evaluate the variance function for the link
 *
 * @param x pointer to a glmer2 object
 * @param var pointer to positions to hold computed values
 *
 * @return var
 */
static double *glmer_var(SEXP x, double *var)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, n = dims[n_POS], fltype = dims[famType_POS];
    double *mu = REAL(GET_SLOT(x, lme4_muSym));

    switch(fltype) {
    case 1: 			/* binomial family with logit or probit link */
    case 2:
	for (i = 0; i < n; i++) {
	    double mui = mu[i];
	    var[i] = mui * (1 - mui);
	}
	break;
    case 3:			/* Poisson with log link */
	for (i = 0; i < n; i++) {
	    var[i] = mu[i];
	}
	break;
    default:
	error(_("General form of glmer_var not yet written"));
    }
    return var;
}

/**
 * Evaluate the derivative of mu wrt eta for the link
 *
 * @param x pointer to a glmer2 object
 * @param dmu_deta pointer to positions to hold computed values
 *
 * @return dmu_deta
 */
static double *glmer_dmu_deta(SEXP x, double *dmu_deta)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, n = dims[n_POS], fltype = dims[famType_POS];
    double *eta = REAL(GET_SLOT(x, lme4_etaSym));

    switch(fltype) {
    case 1: 			/* binomial with logit link */
	for (i = 0; i < n; i++) {
	    double etai = eta[i];
	    double opexp = 1 + exp(etai);
	    
	    dmu_deta[i] = (etai > LTHRESH || etai < MLTHRESH) ?
		DOUBLE_EPS : exp(etai)/(opexp * opexp);
	}
	break;
    case 2:			/* binomial with probit link */
	for (i = 0; i < n; i++) {
	    double tmp = dnorm4(eta[i], 0, 1, 0);
	    dmu_deta[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
	}
	break;
    case 3:			/* Poisson with log link */
	for (i = 0; i < n; i++) {
	    double tmp = exp(eta[i]);
	    dmu_deta[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
	}
	break;
/*     default: { */
/* 	SEXP ans = PROTECT(eval_check(mu_eta, rho, REALSXP, n)); */
/* 	Memcpy(dmu_deta, REAL(ans), n); */
/* 	UNPROTECT(1); */
/*     } */
    }
    return dmu_deta;
}

/* FIXME: Should this update the Sdr_POS value in the deviance slot directly? */
/**
 * Evaluate the deviance residuals
 *
 * @param x pointer to a glmer2 object
 * @param dev_res pointer to an area to hold the result
 *
 * @return sum of the deviance residuals
 */
static double glmer_dev_resids(SEXP x, double *dev_res)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, fltype = dims[famType_POS], n = dims[n_POS];
    double *mu = REAL(GET_SLOT(x, lme4_muSym)),
	*wts = REAL(GET_SLOT(x, install("pwts"))),
	*y = REAL(GET_SLOT(x, lme4_ySym)), sum;

    switch(fltype) {
    case 1: 			/* binomial with logit or probit link */
    case 2:
	for (i = 0; i < n; i++) {
	    double mui = mu[i], yi = y[i];
	    
	    dev_res[i] = 2 * wts[i] *
		(y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
	}
	break;
    case 3:			/* Poisson with log link */
	for (i = 0; i < n; i++) {
	    double mui = mu[i], yi = y[i];
	    dev_res[i] = 2 * wts[i] * (y_log_y(yi, mui) - (yi - mui));
	}
	break;
    default:
	error(_("General form of glmer_dev_resids not yet written"));
    }
    for (i = 0, sum = 0; i < n; i++) sum += dev_res[i];
    return sum;
}

/**
 * Evaluate the linear predictor as model offset + X \beta + Z b
 *
 * @param x pointer to a glmer2 object
 *
 * @return R_NilValue
 */
SEXP glmer_eta(SEXP x)
{
    SEXP moff = GET_SLOT(x, install("moff")),
	fixef = GET_SLOT(x, lme4_fixefSym),
	etap = GET_SLOT(x, lme4_etaSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int ione = 1, n = dims[n_POS], nfe = LENGTH(fixef), q = dims[q_POS];
    /* For this model LENGTH(fixef) != dims[p_POS] == 0 */
    CHM_SP ZXyt = AS_CHM_SP(GET_SLOT(x, lme4_ZXytSym));
    double *eta = REAL(etap), *re = Alloca(ZXyt->nrow, double), one[] = {1,0};
    CHM_DN ceta = AS_CHM_DN(etap), cre = N_AS_CHM_DN(re, ZXyt->nrow);
    R_CheckStack();

    if (LENGTH(moff)) /* initialize eta to model offset */
	Memcpy(eta, REAL(moff), n);
    else
	AZERO(eta, n);
		     /* Add fixed-effects contribution to eta */
    F77_CALL(dgemv)("N", &n, &nfe, one, REAL(GET_SLOT(x, lme4_XSym)), &n,
		    REAL(fixef), &ione, one, eta, &ione);
    lmer2_update_effects(x);	/* evaluate b */
    Memcpy(re, REAL(GET_SLOT(x, lme4_ranefSym)), q);
    re[q] = 0.;
    if (!M_cholmod_sdmult(ZXyt, 1 /* trans */, one, one, cre, ceta, &c))
	error(_("cholmod_sdmult error returned"));
    return R_NilValue;
}

/**
 * Multiply a vector by the virtual T and S matrices represented by ST
 *
 * @param Gp vector of group pointers
 * @param ST compounded matrices S* and T*
 * @param b vector of random effects to be transformed from b* to b
 *
 * @return b after transformation
 */
static double *TS_mult(const int *Gp, SEXP ST, double *b)
{
    int i, ione = 1, nf = LENGTH(ST);

    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int j, k, ncip1 = nci + 1,
	    nlev = (Gp[i + 1] - Gp[i])/nci;

	for (j = 0; j < nlev; j++) {
	    int base = Gp[i] + j * nci;
	    for (k = 0; k < nci; k++) /* multiply by S_i */
		b[base + k] *= st[k * ncip1];
	    if (nci > 1) 	/* multiply by T_i */
		F77_CALL(dtrmv)("L", "N", "U", &nci, st, &nci,
				b + base, &ione);
	}
    }
    return b;
}

/**
 * Evaluate the effects in an lmer2 representation
 *
 * @param L factorization
 *
 * @return cholmod_dense object with \hat{b*} in the first q positions
 * and \hat{\beta} in the next p positions.
 */
static cholmod_dense
*internal_lmer2_effects(cholmod_factor *L)
{
    int i, nt = (L->n);
    cholmod_dense *X,
	*B = M_cholmod_allocate_dense((size_t)nt, 1, (size_t)nt, CHOLMOD_REAL, &c);
    double *bx = (double*)(B->x);
    
    for (i = 0; i < nt; i++) bx[i] = 0;
    if (L->is_super) {
	int ns = (L->nsuper);
	int nr = ((int *)(L->pi))[ns] - ((int *)(L->pi))[ns - 1],
	    nc = ((int *)(L->super))[ns] - ((int *)(L->super))[ns - 1];
	double *x = (double *)(L->x) + ((int *)(L->px))[ns - 1];

	bx[nt - 1] = x[(nc - 1) * (nr + 1)];
    } else {
	bx[nt - 1] = (L->is_ll) ? ((double*)(L->x))[((int*)(L->p))[nt - 1]] : 1;
    }
    if (!(X = M_cholmod_solve(CHOLMOD_Lt, L, B, &c)))
	error(_("cholmod_solve (CHOLMOD_Lt) failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    M_cholmod_free_dense(&B, &c);
    if (!(B = M_cholmod_solve(CHOLMOD_Pt, L, X, &c)))
	error(_("cholmod_solve (CHOLMOD_Pt) failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    M_cholmod_free_dense(&X, &c);
    return B;
}

/**
 * Evaluate the logarithm of the square of the determinant of selected
 * sections of a sparse Cholesky factor.
 *
 * @param ans vector of doubles of sufficient length to hold the result
 * @param nans number of values to calculate
 * @param c vector of length nans+1 containing the cut points
 * @param F factorization
 *
 * @return ans
 */
static double*
chm_log_abs_det2(double *ans, int nans, const int *c, const cholmod_factor *F)
{
    int i, ii  = 0, jj = 0;
    for (i = 0; i < nans; i++) ans[i] = 0;
    if (F->is_super) {
	for (i = 0; i < F->nsuper; i++) {
	    int j, nrp1 = 1 + ((int *)(F->pi))[i + 1] - ((int *)(F->pi))[i],
		nc = ((int *)(F->super))[i + 1] - ((int *)(F->super))[i];
	    double *x = (double *)(F->x) + ((int *)(F->px))[i];

	    for (j = 0; j < nc; j++) {
		int col = j + jj;
		if (col < c[ii]) continue;
		while (col >= c[ii + 1] && ++ii < nans) {};
		if (ii >= nans) break;
		ans[ii] += 2 * log(fabs(x[j * nrp1]));
	    }
	    jj += nc;
	}
    } else {
	int *fi = (int*)(F->i), *fp = (int*)(F->p), j, k;
	double *fx = (double *)(F->x);
	
	for (j = 0; ii < nans && j < F->n; j++) {
	    if (j < c[ii]) continue;
	    for (k = fp[j]; fi[k] != j && k < fp[j + 1]; k++) {};
	    if (fi[k] != j) break; /* what happened to the diagonal element? */
	    while (j >= c[ii + 1] && ++ii < nans) {};
	    if (ii >= nans) break;
	    ans[ii] += log(fx[k] * ((F->is_ll) ? fx[k] : 1.));
	}
    }
    return ans;
}

/**
 * Evaluate the elements of the deviance slot given a factorization of
 * A* and the dimensions vector.
 *
 * @param d pointer to the contents of the slot
 * @param dims dimensions
 * @param L factor of the current A*
 *
 * @return d
 */
static double*
internal_deviance(double *d, const int *dims, const cholmod_factor *L)
{
    int n = dims[n_POS], p = dims[p_POS], q = dims[q_POS];
    int c[] = {0,  q, p + q, p + q + 1};
    double dn = (double) n, dnmp = (double)(n - p);
    
    chm_log_abs_det2(d + ldZ_POS, lr2_POS + 1 - ldZ_POS, c, L);
    d[ML_POS] = d[ldZ_POS] + dn * (1. + d[lr2_POS] + log(2. * PI / dn));
    d[REML_POS] = d[ldZ_POS] + d[ldX_POS] + dnmp *
	(1. + d[lr2_POS] + log(2. * PI / dnmp));
    d[bqd_POS] = d[Sdr_POS] = 0.;
    return d;
}

/**
 * Update A to A* and evaluate its numeric factorization in L.
 *
 * @param deviance Hold the result
 * @param dims dimensions
 * @param Gp length nf+3 vector of group pointers for the rows of A
 * @param ST pointers to the nf ST factorizations of the diagonal
 *     elements of Sigma 
 * @param A symmetric matrix of size Gp[nf+2]
 * @param F factorization to be modified
 *
 */
static void
internal_update_L(double *deviance, int *dims, const int *Gp,
		  SEXP ST, cholmod_sparse *A, cholmod_factor *L)
{
    cholmod_sparse *Ac = M_cholmod_copy_sparse(A, &c);
    int *ai = (int *)(Ac->i), *ap = (int *)(Ac->p), nf = *dims,
	i, j, ione = 1;
    double *ax = (double*)(Ac->x) , one[] = {1, 0};
    

    if ((!Ac->sorted) || Ac->stype <= 0) {
	M_cholmod_free_sparse(&Ac, &c);
	error(_("A must be a sorted cholmod_sparse object with stype > 0"));
    }

    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int base = Gp[i], k, kk;
	int ncip1 = nci + 1, nlev = (Gp[i + 1] - Gp[i])/nci;

	/* if nci == 1 then T == I and the next section is skipped */
	if (nci > 1) {	/* evaluate ith section of SAST */
	    int maxrows = -1;
	    double *db = Calloc(nci * nci, double), /* diagonal blcok */
		*wrk = (double *) NULL;
	    
/* FIXME: calculate and store maxrows in lmer2_create */
	    if (nci > 2) {	/* calculate scratch storage size */
		for (j = 0; j < nlev; j++) {
		    int cj = base + j * nci; /* first column in this group */
		    int nnzm1 = ap[cj + 1] - ap[cj] - 1;
		    if (nnzm1 > maxrows) maxrows = nnzm1;
		}
		wrk = Calloc(maxrows * nci, double);
	    }
	    
	    for (j = 0; j < nlev; j++) {
		int cj = base + j * nci; /* first column in this block */
		int nnz = ap[cj + 1] - ap[cj]; /* nonzeros in column cj */
		int nnzm1 = nnz - 1;
		
		if (nnzm1) { /* elements above diagonal block to update */
		    if (nci == 2)
			F77_CALL(dtrmm)("R", "L", "N", "U", &nnzm1,
					&nci, one, st, &nci,
					ax + ap[cj], &nnz);
		    else {	
			for (k = 0; k < nci; k++) /* copy columns to wrk */
			    Memcpy(wrk + k * maxrows, ax + ap[cj + k],
				   nnzm1);
			F77_CALL(dtrmm)("R", "L", "N", "U", &nnzm1,
					&nci, one, st, &nci, wrk,
					&maxrows);
			for (k = 0; k < nci; k++) /* copy results back */
			    Memcpy(ax + ap[cj + k], wrk + k * maxrows,
				   nnzm1);
		    }
		    /* evaluate T'A for rows and columns in this block */
		    for (k = 0; k < nci; k++) {
			for (kk = 0; kk < nnzm1;) {
			    int ind = ap[cj + k] + kk;
			    if (Gp[i] <= ai[ind]) {
				F77_CALL(dtrmv)("L", "T", "U", &nci,
						st, &nci, ax + ind,
						&ione);
				kk += nci;
			    } else kk++;
			}
		    }
		}
				/* update the diagonal block */
		for (k = 0; k < nci; k++) /* copy upper triangle */
		    Memcpy(db + k * nci, ax + ap[cj + k] + nnzm1, k + 1);
		for (k = 1; k < nci; k++) /* symmetrize */
		    for (kk = 0; kk < k; kk++)
			db[k + kk * nci] = db[kk + k * nci];
		F77_CALL(dtrmm)("L", "L", "T", "U", &nci, &nci, one,
				st, &nci, db, &nci);
		F77_CALL(dtrmm)("R", "L", "N", "U", &nci, &nci, one,
				st, &nci, db, &nci);
		for (k = 0; k < nci; k++) /* restore updated upper triangle */
		    Memcpy(ax + ap[cj + k] + nnzm1, db + k * nci, k + 1);
	    }
	    if (nci > 2) Free(wrk);
				/* evaluate T'AT for all blocks to the right */
	    for (j = Gp[i+1]; j < Gp[nf + 2]; j++) {
		for (k = ap[j]; k < ap[j + 1]; ) {
		    if (ai[k] >= Gp[i + 1]) break;
		    if (ai[k] < Gp[i]) {
			k++;
			continue;
		    }
		    F77_CALL(dtrmv)("L", "T", "U", &nci, st, &nci,
				    ax + k, &ione);
		    k += nci;
		}
	    }
	    Free(db);
	}
				/* Multiply by S from left. */
	for (j = Gp[i]; j < A->ncol; j++)
	    for (k = ap[j]; k < ap[j+1]; k++) {
		int row = ai[k];
		if (row < Gp[i]) continue;
		if (Gp[i + 1] <= row) break;
		ax[k] *= st[((row - Gp[i]) % nci) * ncip1];
	    }
				/* Multiply by S from right */
	for (j = Gp[i]; j < Gp[i + 1]; j += nci) {
	    for (k = 0; k < nci; k++)
		for (kk = ap[j + k]; kk < ap[j + k + 1]; kk++)
		    ax[kk] *= st[k * ncip1];
	}
				/* Increment diagonal */
	for (j = Gp[i]; j < Gp[i + 1]; j++) {
	    k = ap[j + 1] - 1;
	    if (ai[k] != j) error(_("Logic error"));
	    ax[k]++;
	}
    }
				/* force LL' decomp and sorted row indices */
    i = c.final_ll; c.final_ll = TRUE;
    j = c.final_monotonic; c.final_monotonic = TRUE;
    if (!M_cholmod_factorize(Ac, L, &c)) { 
	error(_("cholmod_factorize failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    }	
				/* restore previous settings */
    c.final_ll = i; c.final_monotonic = j;
    internal_deviance(deviance, dims, L);
    M_cholmod_free_sparse(&Ac, &c);
}

/**
 * Return the group in the (nf, Gp) combination to which ind belongs
 *
 * @param ind a row number
 * @param nf number of groups of rows
 * @param Gp group pointers
 *
 */
static int R_INLINE Gp_grp(int ind, int nf, const int *Gp)
{
    int i;
    for (i = 0; i < nf; i++) if (ind < Gp[i + 1]) return i;
    error(_("invalid row index %d (max is %d)"), ind, Gp[nf]);
    return -1;			/* -Wall */
}

/**
 * Update the Vt slot in a glmer or nlmer object.
 *
 * @param x pointer to a glmer or nlmer object
 *
 */
SEXP nlmer_update_Vt(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym),
	Vt = GET_SLOT(x, install("Vt")),
	Zt =  GET_SLOT(x, lme4_ZtSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*vi = INTEGER(GET_SLOT(Vt, lme4_iSym)),
	*vp = INTEGER(GET_SLOT(Vt, lme4_pSym)),
	*zi = INTEGER(GET_SLOT(Zt, lme4_iSym)),
	*zp = INTEGER(GET_SLOT(Zt, lme4_pSym)),
	i, ione = 1, iv, iz, j, mnc, nf = LENGTH(ST),
	ncol = INTEGER(GET_SLOT(Zt, lme4_DimSym))[1];
    int *nc = Calloc(nf, int);
    double **st = Calloc(nf, double*), *tmp,
	*vx = REAL(GET_SLOT(Vt, lme4_xSym)),
	*zx = REAL(GET_SLOT(Zt, lme4_xSym));

    for (i = 0, mnc = 0; i < nf; i++) {	/* populate nc and st */
	SEXP STi = VECTOR_ELT(ST, i);
	nc[i] = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	if (nc[i] > mnc) mnc = nc[i]; /* max num of cols */
	st[i] = REAL(STi);
    }
    tmp = Calloc(mnc, double);
    for (j = 0; j < ncol; j++) {
	int iz2 = zp[j + 1];
				/* premultiply by T' */
	for (iz = zp[j], iv = vp[j]; iz < iz2; iz++) {
	    int k = Gp_grp(zi[iz], nf, Gp);
	    if (nc[k] > 1) {
		int itmp = (zi[iz] - Gp[k]) % nc[k];
		AZERO(tmp, mnc);
		tmp[itmp] = zx[iz];
		for (i = 1; i < nc[k] && (iz + 1) < iz2; i++) {
		    if (zi[iz + 1] != zi[iz] + 1) break;
		    tmp[itmp++] = zx[++iz];
		}
		F77_CALL(dtrmv)("L", "T", "U", &(nc[k]), st[k],
				&(nc[k]), tmp, &ione);
		for (i = 0; i < nc[k] && iv < vp[j + 1]; i++, iv++) {
		    vx[iv] = tmp[i];
		    if (iv == vp[j + 1] - 1) break;
		    if (vi[iv + 1] != vi[iv] + 1) break;
		}
	    } else vx[iv++] = zx[iz++];
	}
	for (iv = vp[j]; iv < vp[j + 1]; iv++) {
	    int k = Gp_grp(vi[iv], nf, Gp);
	    vx[iv] *= st[k][((vi[iv] - Gp[k]) % nc[k]) * (nc[k] + 1)];
	}
    }    
    Free(st); Free(nc); Free(tmp);
    return R_NilValue;
}

/**
 * Update the ranef slot, b=TSu, in a glmer or nlmer object.
 *
 * @param x pointer to a glmer or nlmer object
 *
 */
SEXP nlmer_update_ranef(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), i, ione = 1;
    double *b = REAL(GET_SLOT(x, lme4_ranefSym)),
	*u = REAL(GET_SLOT(x, install("uvec")));
    
    for (i = 0; i < dims[nf_POS]; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *sti = REAL(STi);
	int base = Gp[i], j, k,
	    nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	
	for (j = base; j < Gp[i+1]; j += nci) {
	    for (k = 0; k < nci; k++) { /* premultiply  by S */
		int jj = j + k;
		b[jj] = u[jj] * sti[k];
	    }
	    if (nci > 1)	/* premultiply  by T */
		F77_CALL(dtrmv)("L", "N", "U", &nci, sti,
				&nci, &(u[j]), &ione);
	}
    }
    return R_NilValue;
}

/**
 * Evaluate starting estimates for the elements of ST
 *
 * @param ST pointers to the nf ST factorizations of the diagonal
 *     elements of Sigma 
 * @param Gp length nf+3 vector of group pointers for the rows of A
 * @param A symmetric matrix of size Gp[nf+2]
 *
 */
static void
internal_lmer2_initial(SEXP ST, int *Gp, cholmod_sparse *A)
{
    int *ai = (int*)(A->i), *ap = (int*)(A->p), i, nf = LENGTH(ST);
    double *ax = (double*)(A->x);
    
    if (!(A->sorted) || (A->stype <= 0))
	error(_("A should be upper triangular and sorted"));
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int bb = Gp[i], j, k;
	int ncip1 = nci + 1, nlev = (Gp[i + 1] - bb)/nci;
	
	AZERO(st, nci * nci);
	for (j = 0; j < nlev; j++) {
	    int base = bb + j * nci;
	    for (k = 0; k < nci; k++) {
		int cl = ap[base + k + 1] - 1; /* last element in the column */
		if (ai[cl] != (base + k))
		    error(_("Weird structure in A.  Expecting row %d, got %d"),
			  base + k, ai[cl]);
		st[k * ncip1] += ax[cl];
	    }
	}
	for (k = 0; k < nci; k++)
	    st[k * ncip1] = sqrt(((double) nlev)/(0.375*st[k * ncip1]));
    }
}


static R_INLINE void
make_cholmod_sparse_sorted(cholmod_sparse *A)
{
    if(!A->sorted) {
	int i = M_cholmod_sort(A, &c); 
	if(!i)
	    error(_("cholmod_sort returned error code %d"),i);
    }
}

static void
internal_update_A(cholmod_sparse *ZXyt, SEXP wtP, SEXP offP,
		  cholmod_sparse *A)
{
    cholmod_sparse *ts1 = M_cholmod_copy_sparse(ZXyt, &c), *ts2;
    int *ai, *ap, *zi, *zp, i, j, m = ts1->nrow, n = ts1->ncol,
	ol = LENGTH(offP), wl = LENGTH(wtP);
    double *ax, *zx, *off = REAL(offP), *wts = REAL(wtP);

    make_cholmod_sparse_sorted(ts1);
    zi = (int*)(ts1->i); zp = (int*)(ts1->p);
    zx  = (double*)(ts1->x);
    if (ol) {		/* skip if length 0 */
	if (ol != n) {
	    M_cholmod_free_sparse(&ts1, &c);
	    error(_("Length of offset is %d, should be %d"), ol, n);
	}
	for (j = 0; j < n; j++) { /* iterate over columns */
	    int ind = zp[j + 1] - 1;
	    if (zi[ind] != (m - 1)) {
		M_cholmod_free_sparse(&ts1, &c);
		error(_("missing y position in ZXyt at column %d"), j+1);
	    }
	    zx[ind] += off[j];	/* *add* offset to -y */
	}
    }
    if (wl) {		/* skip if length 0 */
	if (wl != n) {
	    M_cholmod_free_sparse(&ts1, &c);
	    error(_("Length of offset is %d, should be %d"), wl, n);
	}
	for (j = 0; j < n; j++) { /* iterate over columns */
	    double wt = sqrt(wts[j]);
	    for (i = zp[j]; i < zp[j + 1]; i++) zx[i] *= wt;
	}
    }
    
    ts2 = M_cholmod_aat(ts1, (int*)NULL, (size_t)0, 1, &c);
    M_cholmod_free_sparse(&ts1, &c);
    ts1 = M_cholmod_copy(ts2, +1/*upper triangle*/, +1/*values*/, &c);
    M_cholmod_free_sparse(&ts2, &c);
    make_cholmod_sparse_sorted(ts1);
    make_cholmod_sparse_sorted(A);
    ai = (int*)(A->i); zi = (int*)(ts1->i);
    ap = (int*)(A->p); zp = (int*)(ts1->p);
    ax = (double*)(A->x); zx = (double*)(ts1->x);
    for (j = 0; j < m; j++) {
	if (ap[j + 1] != zp[j + 1]) {
	    M_cholmod_free_sparse(&ts1, &c);
	    error(_("pattern mismatch for A and update at column %d"),
		  j + 1);
	}
	for (i = ap[j]; i < ap[j + 1]; i++) {
	    if (ai[i] != zi[i]) {
		M_cholmod_free_sparse(&ts1, &c);
		error(_("pattern mismatch for A and update at column %d"),
		      j + 1);
	    }
	    ax[i] = zx[i];
	}
    }
    M_cholmod_free_sparse(&ts1, &c);
}

/**
 * Evaluate new weights and working residuals.
 *
 * @param x pointer to a glmer2 object
 */
SEXP glmer_reweight(SEXP x)
{
    SEXP fixef = GET_SLOT(x, lme4_fixefSym),
	off = GET_SLOT(x, lme4_offsetSym),
	moff = GET_SLOT(x, install("moff")),
	wts = GET_SLOT(x, lme4_weightsSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, ione = 1, lo = LENGTH(moff),
	n = dims[n_POS], nfe = LENGTH(fixef);
    double
	*dmu_deta = Calloc(n, double),
	*eta = REAL(GET_SLOT(x, lme4_etaSym)),
	*mo = REAL(moff), *mu = REAL(GET_SLOT(x, lme4_muSym)),
	*var = Calloc(n, double), *w = REAL(wts),
	*y = REAL(GET_SLOT(x, lme4_ySym)), *z = REAL(off),
	one[] = {1, 0};
    CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym)),
	ZXyt = AS_CHM_SP(GET_SLOT(x, lme4_ZXytSym));
    R_CheckStack();      
				/* initialize weights to prior wts */
    Memcpy(w, REAL(GET_SLOT(x, install("pwts"))), n); 
    glmer_linkinv(x);	      /* evaluate mu */
    glmer_dmu_deta(x, dmu_deta);
    glmer_var(x, var);
    for (i = 0; i < n; i++) {
	w[i] *= dmu_deta[i] * dmu_deta[i]/var[i];
	/* store negative of adj. wrk. variate in offset */
 	z[i] = (lo ? mo[i] : 0) - eta[i] - (y[i] - mu[i])/dmu_deta[i];
    }
    /* Add the contribution of X\beta to the offset. This is a kludge. */
    F77_CALL(dgemv)("N", &n, &nfe, one, REAL(GET_SLOT(x, lme4_XSym)),
		    &n, REAL(fixef), &ione, one, z, &ione);
    internal_update_A(ZXyt, wts, off, A);
    Free(dmu_deta); Free(var);
    return R_NilValue;
}

/**
 * Extract the parameters from ST
 *
 * @param ST ST slot from an lmer2 object
 * @param pars double vector of the appropriate length
 *
 * @return pointer to the parameter vector
 */
static double
*internal_lmer2_getPars(SEXP ST, double *pars)
{
    int i, nf = LENGTH(ST), pos = 0;
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int j, k, ncp1 = nci + 1;

	for (j = 0; j < nci; j++) pars[pos++] = st[j * ncp1];
	for (j = 0; j < (nci - 1); j++)
	    for (k = j + 1; k < nci; k++)
		pars[pos++] = st[k + j * nci];
    }
    return pars;
}

/**
 * Extract the parameters from the ST slot of an lmer2 object
 *
 * @param x an lmer2 object
 *
 * @return pointer to a REAL vector
 */
SEXP lmer2_getPars(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    SEXP ans = PROTECT(allocVector(REALSXP,
				   INTEGER(GET_SLOT(x, lme4_dimsSym))[np_POS]));

    internal_lmer2_getPars(ST, REAL(ans));
    UNPROTECT(1); 
    return ans;
}

/**
 * Update the ST slot of an lmer2 object
 *
 * @param pars double vector of the appropriate length
 * @param ST ST slot from an lmer2 object
 *
 * @return pointer to the updated ST object
 */
static SEXP
internal_lmer2_setPars(const double *pars, SEXP ST)
{
    int i, nf = LENGTH(ST), pos = 0;
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int j, k, ncp1 = nci + 1;

	for (j = 0; j < nci; j++) st[j * ncp1] = pars[pos++];
	for (j = 0; j < (nci - 1); j++)
	    for (k = j + 1; k < nci; k++)
		st[k + j * nci] = pars[pos++];
    }
    return ST;
}

/**
 * Update the ST slot of an lmer2 object from a REAL vector of
 * parameters and update the Cholesky factorization
 *
 * @param x an lmer2 object
 * @param pars a REAL vector of the appropriate length
 *
 * @return x
 */
SEXP lmer2_setPars(SEXP x, SEXP pars)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int npar = INTEGER(GET_SLOT(x, lme4_dimsSym))[np_POS];
    CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym));
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    R_CheckStack();

    if (!isReal(pars) || LENGTH(pars) != npar)
	error(_("pars must be a real vector of length %d"), npar);
    internal_lmer2_setPars(REAL(pars), ST);
    internal_update_L(REAL(GET_SLOT(x, lme4_devianceSym)),
		      INTEGER(GET_SLOT(x, lme4_dimsSym)), 
		      INTEGER(GET_SLOT(x, lme4_GpSym)), ST, A, L);
    return x;
}

/**
 * Extract the deviance from an lmer2 object
 *
 * @param x an lmer2 object
 * @param which scalar integer < 0 => REML, 0 => native, > 0 => ML
 *
 * @return scalar REAL value
 */
SEXP lmer2_deviance(SEXP x, SEXP which)
{
    int w = asInteger(which);
    int POS = (w < 0 || (!w && isREML(x))) ? REML_POS : ML_POS; 

    return ScalarReal(REAL(GET_SLOT(x, lme4_devianceSym))[POS]);
}

    
/**
 * Update the contents of the fixef and ranef slots
 *
 * @param x an lmer2 object
 *
 * @return R_NilValue
 */
SEXP lmer2_update_effects(SEXP x)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), i;
    int q = dims[q_POS];
    double *b = REAL(GET_SLOT(x, lme4_ranefSym)), *bstbx,
	*dev = REAL(GET_SLOT(x, lme4_devianceSym));
    CHM_DN bstarb;
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    R_CheckStack();

    bstarb = internal_lmer2_effects(L);
    bstbx = (double*)(bstarb->x);
    Memcpy(b, bstbx, q);
    for (i = 0, dev[bqd_POS] = 0; i < q; i++) /* accumulate ssqs of bstar */
	dev[bqd_POS] += bstbx[i] * bstbx[i];
    /* FIXME: apply the permutation when copying */
    Memcpy(REAL(GET_SLOT(x, lme4_fixefSym)), bstbx + q, dims[p_POS]);
    M_cholmod_free_dense(&bstarb, &c);
    TS_mult(INTEGER(GET_SLOT(x, lme4_GpSym)),
	    GET_SLOT(x, lme4_STSym), b);
    return R_NilValue;
}

/**
 * Return the REML or ML conditional estimate of sigma, the standard
 * deviation of the per-observation noise term.
 *
 * @param REML non-zero for REML estimate, 0 for ML estimate
 * @param dims vector of dimensions
 * @param deviance vector of deviance components
 */
static R_INLINE double
internal_lmer2_sigma(int REML, const int* dims, const double* deviance)
{
    return sqrt(exp(deviance[lr2_POS])/
		((double)(dims[n_POS] - (REML ? dims[p_POS] : 0))));
}


/**
 * Extract the estimate of the scale factor from an lmer2 object
 *
 * @param x an lmer2 object
 * @param which scalar integer (< 0 => REML, 0 => native, > 0 => ML)
 *
 * @return scalar REAL value
 */
SEXP lmer2_sigma(SEXP x, SEXP which)
{
    int w = asInteger(which);
		
    return ScalarReal(internal_lmer2_sigma(w < 0 || (!w && isREML(x)),
					  INTEGER(GET_SLOT(x, lme4_dimsSym)),
					  REAL(GET_SLOT(x, lme4_devianceSym))));
}

/**
 * Extract the unscaled lower Cholesky factor of the relative
 * variance-covariance matrix for the fixed-effects in an lmer2 object.
 *
 * @param x an lmer2 object
 *
 * @return a REAL p by p lower triangular matrix (it's a matrix, not a Matrix)
 */

SEXP lmer2_vcov(SEXP x)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), *select;
    int i, p = dims[p_POS], q = dims[q_POS];
    SEXP ans = PROTECT(allocMatrix(REALSXP, p, p));

    if (p) {
	CHM_SP Lsp, Lred; /* sparse and reduced-size sparse */
	CHM_DN Ld;
	CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym)), Lcp;
	select = Alloca(p, int);
	R_CheckStack();

	 /* need a copy of L because factor_to_sparse changes 1st arg */
	Lcp = M_cholmod_copy_factor(L, &c);
	if (!(Lcp->is_ll))
	    if (!M_cholmod_change_factor(Lcp->xtype, 1, 0, 1, 1, Lcp, &c))
		error(_("cholmod_change_factor failed with status %d"), c.status);
	Lsp = M_cholmod_factor_to_sparse(Lcp, &c);
	M_cholmod_free_factor(&Lcp, &c);
	for (i = 0; i < p; i++) select[i] = q + i;
	Lred = M_cholmod_submatrix(Lsp, select, p, select, p,
				 1 /* values */, 1 /* sorted */, &c);
	M_cholmod_free_sparse(&Lsp, &c);
	Ld = M_cholmod_sparse_to_dense(Lred, &c);
	M_cholmod_free_sparse(&Lred, &c);
	Memcpy(REAL(ans), (double*)(Ld->x), p * p);
	M_cholmod_free_dense(&Ld, &c);
/* FIXME: This does not allow for a possible P_X permutation  */
	F77_CALL(dtrtri)("L", "N", &p, REAL(ans), &p, &i);
	if (i)
	    error(_("Lapack routine dtrtri returned error code %d"), i);
    }
    UNPROTECT(1);
    return ans;
}

/**
 * Extract the conditional modes of the random effects as a list of matrices
 *
 * @param x Pointer to an mer object
 *
 * @return a list of matrices containing the conditional modes of the
 * random effects
 */
SEXP lmer2_ranef(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym),
        cnames = GET_SLOT(x, lme4_cnamesSym),
	flist = GET_SLOT(x, lme4_flistSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	i, ii, jj, nf = LENGTH(flist);
    SEXP val = PROTECT(allocVector(VECSXP, nf));
    double *b = REAL(GET_SLOT(x, lme4_ranefSym));

    lmer2_update_effects(x);
    setAttrib(val, R_NamesSymbol,
	      duplicate(getAttrib(flist, R_NamesSymbol)));
    for (i = 0; i < nf; i++) {
	SEXP nms, rnms = getAttrib(VECTOR_ELT(flist, i), R_LevelsSymbol);
	int nci = INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol))[0];
	int mi = length(rnms);
	double *bi = b + Gp[i], *mm;

	SET_VECTOR_ELT(val, i, allocMatrix(REALSXP, mi, nci));
	setAttrib(VECTOR_ELT(val, i), R_DimNamesSymbol, allocVector(VECSXP, 2));
	nms = getAttrib(VECTOR_ELT(val, i), R_DimNamesSymbol);
	SET_VECTOR_ELT(nms, 0, duplicate(rnms));
	SET_VECTOR_ELT(nms, 1, duplicate(VECTOR_ELT(cnames, i)));
	mm = REAL(VECTOR_ELT(val, i));
	for (jj = 0; jj < nci; jj++)
	    for(ii = 0; ii < mi; ii++)
		mm[ii + jj * mi] = bi[jj + ii * nci];
    }
    UNPROTECT(1);
    return val;
}

/**
 * Extract the posterior variances of the random effects
 *
 * @param x pointer to a mer object
 *
 * @return pointer to a list of arrays
 */
SEXP lmer2_postVar(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    double *deviance = REAL(GET_SLOT(x, lme4_devianceSym)), one = 1;
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, j, nf = dims[nf_POS], p = dims[p_POS], q = dims[q_POS];
    int ppq = p + q;
    SEXP ans = PROTECT(allocVector(VECSXP, nf));
    double sc = internal_lmer2_sigma(isREML(x), dims, deviance);
    CHM_SP rhs, B, Bt, BtB;
    CHM_DN BtBd;
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym)), Lcp = (CHM_FR)NULL;
    int *Perm = (int*)(L->Perm), *iperm = Alloca(ppq, int),
	*fset = Alloca(ppq, int);
    R_CheckStack();
    
    for (j = 0; j < ppq; j++) {
	iperm[Perm[j]] = j;
	fset[j] = j;
    }
    if (!L->is_ll) {
	Lcp = M_cholmod_copy_factor(L, &c);
	L = Lcp;
	j = M_cholmod_change_factor(CHOLMOD_REAL, TRUE/*ll*/,
				    FALSE/*super*/, TRUE/*packed*/,
				    TRUE/*sorted*/, L, &c);
	if (!j) error(_("cholmod_change_factor failed"));
    }
    sc = sc * sc;		/* variance scale factor */
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	int j, k, kk, nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int nlev = (Gp[i + 1] - Gp[i])/nci;
	SEXP ansi = PROTECT(alloc3Darray(REALSXP, nci, nci, nlev));
	int ncip1 = nci + 1, ncisqr = nci * nci;
	double *vv = REAL(ansi),
	    *st = Memcpy(Calloc(ncisqr, double), REAL(STi), ncisqr);

	SET_VECTOR_ELT(ans, i, ansi); UNPROTECT(1);
	AZERO(vv, ncisqr * nlev);
	rhs = M_cholmod_allocate_sparse((size_t)(ppq + 1),
					(size_t) nci, (size_t) nci,
					1/*sorted*/, 1/*packed*/,
					0/*stype*/, CHOLMOD_REAL, &c);
	((int*)(rhs->p))[0] = 0;
	for (k = 0; k < nci; k++) {
	    ((double*)(rhs->x))[k] = 1.;
	    ((int*)(rhs->p))[k + 1] = k + 1;
	}
	for (k = 0; k < nci; k++) {
	    double mult = st[k * ncip1];
	    for (kk = k + 1; kk < nci; kk++)
		st[kk + k * nci] *= mult;
	}
	for (j = 0; j < nlev; j++) {
	    int *ip, *pp, base = Gp[i] + j * nci;
	    double *xp;
	    
	    for (k = 0; k < nci; k++)
		((int*)(rhs->i))[k] = iperm[base + k];
	    B = M_cholmod_spsolve(CHOLMOD_L, L, rhs, &c);
	    ip = (int*)(B->i);
	    pp = (int*)(B->p);
	    xp = (double*)(B->x);
	    if (nci == 1) {
		for (k = 0; k < pp[1]; k++)
		    if (ip[k] < ppq) vv[j] += xp[k] * xp[k];
		vv[j] *= sc * st[0] * st[0];
	    } else {
		double *vvj = vv + j * ncisqr;
		Bt = M_cholmod_transpose(B, TRUE/*values*/, &c);
		BtB = M_cholmod_aat(Bt, fset, (size_t)ppq, 1/*mode*/,&c);
		M_cholmod_free_sparse(&Bt, &c);
		BtBd = M_cholmod_sparse_to_dense(BtB, &c);
		M_cholmod_free_sparse(&BtB, &c);
		Memcpy(vvj, (double*)(BtBd->x), ncisqr);
		M_cholmod_free_dense(&BtBd, &c);
		F77_CALL(dtrmm)("L", "L", "N", "N", &nci, &nci,
				&one, st, &nci, vvj, &nci);
		F77_CALL(dtrmm)("R", "L", "T", "N", &nci, &nci,
				&sc, st, &nci, vvj, &nci);
	    }
	    M_cholmod_free_sparse(&B, &c);
	}
	M_cholmod_free_sparse(&rhs, &c);
	Free(st);
    }
    if (L == Lcp) M_cholmod_free_factor(&L, &c);
    UNPROTECT(1);
    return ans;
}

/**
 * Update the coefficients beta and bstar given the current L
 *
 * @param sigma current standard deviation of the noise term
 * @param L current factorization of A*
 * @param bbsthat current conditional modes
 * @param bbstnew array to receive updated coefficients
 *
 * @return squared length of spherical Gaussian deviate
 */
static double
internal_betabst_update(double sigma, cholmod_factor *L,
		       cholmod_dense *bbsthat,
		       cholmod_dense *bbstnew)
{
    int i, ppq1 = L->n;
    int ppq = ppq1 - 1;
    double *bx, *hx = (double*)(bbsthat->x),
	*nx = (double*)(bbstnew->x), ans = 0;
    cholmod_dense *chb;

    nx[ppq] = 0.;
    for (i = 0; i < ppq; i++) {
	double nr = norm_rand();
	ans += nr * nr;
	nx[i] = sigma * nr;
    }
    				/* chb := L^{-T} %*% bnew */
    chb = M_cholmod_solve(CHOLMOD_Lt, L, bbstnew, &c);
    bx = (double*)(chb->x);
    for (i = 0; i < ppq; i++) nx[i] = bx[i] + hx[i];
    M_cholmod_free_dense(&chb, &c);

    return ans;
}

/**
 * Update the relative variance-covariance matrices by sampling from a
 * Wishart distribution with scale factor determined by the current
 * sample of random effects.
 *
 * @param sigma current value of sigma
 * @param Gp vector of group pointers
 * @param bnew current sample from the random effects
 * @param ST list of factorizations to update
 */
static void
internal_ST_update(double sigma, int trans, const int *Gp,
		   const double *bnew, SEXP ST, double *dest)
{
    int i, j, k, info, nf = LENGTH(ST), pos = 0;
    double one = 1, zero = 0, sigsq = sigma * sigma;
    double sigsqinv = 1/sigsq;

    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi), sd;
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int nlev = (Gp[i + 1] - Gp[i])/nci;

	if (nci == 1) {		/* fast update for scalar */
	    double ssq = 0;
	    for (j = 0; j < nlev; j++) ssq += bnew[Gp[i] + j] * bnew[Gp[i] + j];
	    sd = sqrt(ssq/rchisq((double)nlev)) * st[0];
	    st[0] = sd/sigma;
	    dest[pos++] = (trans ? 2 * log(sd) : sd * sd);
	} else {
	    int ncip1 = nci + 1, ncisqr = nci * nci;
	    double *scal = Calloc(ncisqr, double), /* factor of scale matrix */
		*wfac = Calloc(ncisqr, double); /* factor of Wishart variate */
				/* create the scale matrix (close to identity) */
	    AZERO(scal, ncisqr);
	    F77_CALL(dsyrk)("L", "N", &nci, &nlev, &sigsqinv, bnew + Gp[i], &nci,
			    &zero, scal, &nci);
	    F77_CALL(dpotrf)("L", &nci, scal, &nci, &info);
	    if (info) error(_("Crossprod of b*[[%d]] not positive definite"), i + 1);
	    for (i = 0; i < nci; i++) /* premultiply by S */
		for (j = 0; j <= i; j++) scal[i + j * nci] *= st[i * ncip1];
	    F77_CALL(dtrmm)("L", "L", "N", "U", /* premultiply by T */
			    &nci, &nci, &one, st, &nci, scal, &nci);
				/* generate (lower) random factor from std Wishart */
	    std_rWishart_factor((double)(nlev - nci + 1), nci, 0, wfac);
	    F77_CALL(dtrsm)("L", "L", "N", "N",/* form a scaled variance factor */
			    &nci, &nci, &one, wfac, &nci, scal, &nci);
	    Memcpy(st, scal, ncisqr);
	    for (j = 0; j < nci; j++) { /* form the ST representation */
		for (i = j + 1; i < nci ; i++)
		    st[i + j * nci] /= st[j * ncip1];
	    }
				/* Overwrite wfac with variance-covariance */
	    F77_CALL(dsyrk)("L", "N", &nci, &nci, &sigsq, scal, &nci,
			    &zero, wfac, &nci);
	    if (trans) {
		for (j = 0; j < nci; j++) {
		    double vj = wfac[j * ncip1]; /* jth variance */
		    double sdj = sqrt(vj);

		    for (k = 0; k < j; k++) /* jth row in lower tri */
			wfac[k * nci + j] = atanh(wfac[k * nci + j]/sdj);
		    for (k = j + 1; k < nci; k++)
			wfac[j * nci + k] /= sdj; /* jth col in lower tri */
		    wfac[j * ncip1] = log(vj);
		}
	    }
	    for (j = 0; j < nci; j++) dest[pos++] = wfac[j * ncip1];
	    for (j = 1; j < nci; j++) {
		for (k = 0; k < j; k++)
		    dest[pos++] = wfac[k * nci + j];
	    }
	    Free(scal); Free(wfac);
	}
    }
}

/**
 * Generate a Markov-Chain Monte Carlo sample from a fitted
 * linear mixed model.
 *
 * @param x pointer to an lmer2 object
 * @param savebp pointer to a logical scalar indicating if the
 * random-effects should be saved
 * @param nsampp pointer to an integer scalar of the number of samples
 * to generate
 * @param transp pointer to an logical scalar indicating if the
 * variance components should be transformed.
 *
 * @return a matrix
 */
SEXP lmer2_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp,
		  SEXP verbosep, SEXP deviancep)
{
    SEXP ST = GET_SLOT(x, lme4_STSym), Pars = PROTECT(lmer2_getPars(x)), ans;
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)),
	*Gp = INTEGER(GET_SLOT(x, lme4_GpSym));
    int dPOS = dims[REML_POS] ? REML_POS : ML_POS, i, j,
	n = dims[n_POS], nsamp = asInteger(nsampp),
	p = dims[p_POS], q = dims[q_POS],
	saveb = asLogical(savebp), trans = asLogical(transp),
	verbose = asLogical(verbosep), dev = asLogical(deviancep);
    int qpp1 = q + p + 1;
    double
	*deviance = REAL(GET_SLOT(x, lme4_devianceSym)),
	*ansp, df = n - (dims[REML_POS] ? p : 0);
    int nrbase = p + 1 + dims[np_POS]; /* rows always included */
    int nrtot = nrbase + dev + (saveb ? q : 0);
    cholmod_dense *chhat,
	*chnew = M_cholmod_allocate_dense(qpp1, 1, qpp1, CHOLMOD_REAL, &c);
    double *bstar = (double*) NULL, *nx = (double*)(chnew->x);
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym));
    R_CheckStack();

    if (nsamp <= 0) nsamp = 1;
    ans = PROTECT(allocMatrix(REALSXP, nrtot, nsamp));
    ansp = REAL(ans);
    for (i = 0; i < nrtot * nsamp; i++) ansp[i] = NA_REAL;
    if (saveb) bstar = Calloc(q, double);
    GetRNGstate();
    if (verbose) Rprintf("%12s %14s\n", "sigma", "deviance");

    for (i = 0; i < nsamp; i++) {
	double *col = ansp + i * nrtot, sigma;
				/* simulate and store new value of sigma */
	sigma = exp(deviance[lr2_POS]/2)/sqrt(rchisq(df));
	col[p] = (trans ? 2 * log(sigma) : sigma * sigma);
	/* simulate new fixed and random effects */
				/* Evaluate conditional estimates */
	chhat = internal_lmer2_effects(L);
	internal_betabst_update(sigma, L, chhat, chnew);
	M_cholmod_free_dense(&chhat, &c);
				/* Store beta */
	for (j = 0; j < p; j++) col[j] = nx[q + j];
	if (saveb) {		/* Optionally store b */
	    TS_mult(Gp, ST, Memcpy(bstar, nx, q));
	    for (j = 0; j < q; j++) col[nrbase + dev + j] = bstar[j];
	}
	internal_ST_update(sigma, trans, Gp, nx, ST, col + p + 1);
				/* Refactor and evaluate deviance */
	internal_update_L(deviance, dims, Gp, ST, A, L);
				/* store new variance-covariance parameters */
				/* optionally store deviance */
	/* FIXME: This deviance should be for sigma, beta, b, ST */
	if (dev) col[nrbase] = deviance[dPOS]; 
	if (verbose) Rprintf("%12.6g %14.8g\n", sigma, deviance[dPOS]);
    }
    PutRNGstate();
    if (saveb) Free(bstar);
    M_cholmod_free_dense(&chnew, &c);
				/* Restore pars, refactor, etc. */
    lmer2_setPars(x, Pars);
    UNPROTECT(2);
    return ans;
}

SEXP nlmer_validate(SEXP x)
{
    SEXP GpP = GET_SLOT(x, lme4_GpSym),
	ST = GET_SLOT(x, lme4_STSym),
 	devianceP = GET_SLOT(x, lme4_devianceSym),
	dimsP = GET_SLOT(x, lme4_dimsSym),
	fixefP = GET_SLOT(x, lme4_fixefSym),
	flistP = GET_SLOT(x, lme4_flistSym),
	ranefP = GET_SLOT(x, lme4_ranefSym),
	weightsP = GET_SLOT(x, lme4_weightsSym) ;
    int *Gp = INTEGER(GpP), *dd = INTEGER(dimsP);
    int nf = dd[nf_POS], n = dd[n_POS], p = dd[p_POS], q = dd[q_POS], s = dd[s_POS];
    CHM_SP Xt = AS_CHM_SP(GET_SLOT(x, install("Xt"))),
	Zt = AS_CHM_SP(GET_SLOT(x, install("Zt")));
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    R_CheckStack();

    if (!LENGTH(devianceP))
	return mkString(_("deviance slot must have positive length"));
    if (nf < 1 || LENGTH(flistP) != nf || LENGTH(ST) != nf)
	return mkString(_("Slots ST, and flist must have length nf"));
    if (LENGTH(GpP) != (nf + 1))
	return mkString(_("Slot Gp must have length nf + 1"));
    if (Gp[0] != 0 || Gp[nf] != q)
	return mkString(_("Gp[1] != 0 or Gp[nf+1] != q"));
    if (LENGTH(ranefP) != q)
	return mkString(_("Slot ranef must have length q"));
    if (LENGTH(fixefP) != p)
	return mkString(_("Slot fixef must have length p"));
    if (LENGTH(weightsP) && LENGTH(weightsP) != n)
	return mkString(_("Slot weights must have length 0 or n"));
    if (Zt->nrow != q || Zt->ncol != n * s)
	return mkString(_("Slot Zt must have dimensions q by n*s"));
    if (Xt->nrow != p || Xt->ncol != n * s)
	return mkString(_("Slot Xt must have dimensions p by n*s"));
    if (L->n != q || !L->is_ll || !L->is_monotonic)
	return mkString(_("Slot L must be a monotonic LL' factorization of size q"));
    return ScalarLogical(1);
}

/**
 * Check validity of an lmer2 object.
 *
 * @param x Pointer to an lmer2 object
 *
 * @return TRUE if the object is a valid lmer2 object, else a string
 * describing the nature of the violation.
 */
SEXP lmer2_validate(SEXP x)
{
    SEXP GpP = GET_SLOT(x, lme4_GpSym),
	ST = GET_SLOT(x, lme4_STSym),
	devianceP = GET_SLOT(x, lme4_devianceSym),
	dimsP = GET_SLOT(x, lme4_dimsSym),
	fixefP = GET_SLOT(x, lme4_fixefSym),
	flistP = GET_SLOT(x, lme4_flistSym),
	offsetP = GET_SLOT(x, lme4_offsetSym),
	ranefP = GET_SLOT(x, lme4_ranefSym),
	weightsP = GET_SLOT(x, lme4_weightsSym) ;
    int *Gp = INTEGER(GpP), *dd = INTEGER(dimsP);
    int nf = dd[nf_POS], n = dd[n_POS], p = dd[p_POS], q = dd[q_POS];
    int i, nq, ppq1 = p + q + 1;
    CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym)),
	ZXyt = AS_CHM_SP(GET_SLOT(x, lme4_ZXytSym));
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    R_CheckStack();
				/* check lengths */
    if (nf < 1 || LENGTH(flistP) != nf || LENGTH(ST) != nf)
	return mkString(_("Slots ST, and flist must have length nf"));
    if (LENGTH(GpP) != (nf + 3))
	return mkString(_("Slot Gp must have length nf + 3"));
    if (Gp[0] != 0 || Gp[nf + 2] != ppq1)
	return mkString(_("Gp[1] != 0 or Gp[nf+3] != p+q+1"));
    if (p && LENGTH(fixefP) != p) /* p == 0 is a special case */
	return mkString(_("Slot fixef must have length p"));
    if (LENGTH(ranefP) != q)
	return mkString(_("Slot ranef must have length q"));
    if (LENGTH(weightsP) && LENGTH(weightsP) != n)
	return mkString(_("Slot weights must have length 0 or n"));
    if (LENGTH(offsetP) && LENGTH(offsetP) != n)
	return mkString(_("Slot offset must have length 0 or n"));
    if (LENGTH(devianceP) != (Sdr_POS + 1) ||
	LENGTH(getAttrib(devianceP, R_NamesSymbol)) != (Sdr_POS + 1))
	return mkString(_("deviance slot not named or incorrect length"));
    if (ZXyt->nrow != ppq1 || ZXyt->ncol != n)
	return mkString(_("Slot ZXyt must have dimensions (p+q+1) by n"));
    if (A->nrow != ppq1 || A->ncol != ppq1 || A->stype <= 0)
	return mkString(_("Slot A must be symmetric (upper) of size (p+q+1)"));
    if (L->n != ppq1 || !L->is_ll || !L->is_monotonic)
	return mkString(_("Slot L must be a monotonic LL' factorization of size (p+q+1)"));

    nq = 0;
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i), fli = VECTOR_ELT(flistP, i);
	int *dm = INTEGER(getAttrib(STi, R_DimSymbol));
	if (!isMatrix(STi) || !isReal(STi) || dm[0] != dm[1])
	    return
		mkString(_("Slot ST must be a list of square numeric matrices"));
	if (Gp[i] > Gp[i + 1])
	    return mkString(_("Gp must be non-decreasing"));
	if (!isFactor(fli))
	    return mkString(_("flist must be a list of factors"));
	nq += dm[0] * LENGTH(getAttrib(fli, R_LevelsSymbol));
    }
    if (q != nq)
	return mkString(_("q is not sum of columns by levels"));

    return ScalarLogical(1);
}

static int flType(SEXP family)
{
    const char *fam = CHAR(asChar(getListElement(family, "family"))),
	*lnk = CHAR(asChar(getListElement(family, "link")));

    if ((!strcmp(fam, "gaussian")) && (!strcmp(lnk, "identity")))
	return -1;
    if ((!strcmp(fam, "binomial")) && (!strcmp(lnk, "logit")))
	return 1;
    if ((!strcmp(fam, "binomial")) && (!strcmp(lnk, "probit")))
	return 2;
    if ((!strcmp(fam, "poisson")) && (!strcmp(lnk, "log")))
	return 3;
    return 0;
}

static SEXP NullFrame(void)
{
    SEXP ans = PROTECT(allocVector(VECSXP, 0)); /* empty list */
    setAttrib(ans, R_ClassSymbol, mkString("data.frame"));
    UNPROTECT(1);
    return ans;
}

/**
 * Check Zt for nesting and create the initial permutation
 *
 * @param dims
 * @param Zt
 *
 * @return x without the names attribute
 */
static void Zt_perm(const int *Gp, cholmod_sparse *Zt, SEXP ST,
		    int *dims, int *Perm)
{
    int i;
    cholmod_sparse *ts1, *ts2;
    cholmod_factor *L;

    for (i = 0; i < dims[q_POS]; i++) Perm[i] = i;
    dims[isNest_POS] = TRUE;
    if (dims[nf_POS] > 1) {
	ts1 = M_cholmod_aat(Zt, (int*)NULL/* fset */,(size_t)0,
			    CHOLMOD_PATTERN, &c);
	ts2 = M_cholmod_copy(ts1, -1/*lower triangle*/, CHOLMOD_PATTERN, &c);
	M_cholmod_free_sparse(&ts1, &c);
	if (!check_nesting(dims[nf_POS], ST, Gp, (int*)(ts2->p))) {
	    dims[isNest_POS] = FALSE;
	    L = M_cholmod_analyze(Zt, &c);
	    if (!L)
		error(_("cholmod_analyze returned with status %d"), c.status);
	    Memcpy(Perm, (int*)(L->Perm), dims[q_POS]);
	    M_cholmod_free_factor(&L, &c);
	}
	M_cholmod_free_sparse(&ts2, &c);
    }
}

static void ZXyt_create(SEXP Ztl, SEXP Xp, SEXP yp, SEXP val)
{
    int *Gp = INTEGER(GET_SLOT(val, lme4_GpSym)), *Perm, 
	*dims = INTEGER(GET_SLOT(val, lme4_dimsSym)),
	*i1, *i2, *p1, *p2, i, j;
    int hasY = LENGTH(yp), n = dims[n_POS],
	nf = dims[nf_POS], p = dims[p_POS], q = dims[q_POS];
    double *dest, *src, *y = REAL(yp);
    CHM_DN Xy;
    CHM_FR L;
    SEXP ST = GET_SLOT(val, lme4_STSym);
    CHM_SP A, Zt = AS_CHM_SP(VECTOR_ELT(Ztl, 0)), ts1, ts2, ts3;
    R_CheckStack();

    if (Zt->ncol != n)
	error(_("dimension mismatch: ncol(Ztl[[1]]) = %d != length(y) = %d"),
	      Zt->ncol, n);
    if (Zt->nrow != Gp[1])
	error(_("Expected nrow(Ztl[[1]]) to be %d, got %d"), Gp[1], Zt->nrow);
    ts1 = M_cholmod_copy_sparse(Zt, &c); /* use cholmod storage, not R */
    Zt = ts1;
    for (i = 1; i < nf; i++) {
	ts1 = AS_CHM_SP(VECTOR_ELT(Ztl, i));
	R_CheckStack();
	if (ts1->ncol != n)
	    error(_("dimension mismatch: ncol(Ztl[[%d]]) = %d != length(y) = %d"),
		  i + 1, ts1->ncol, n);
	if (ts1->nrow != (Gp[i + 1] - Gp[i]))
	    error(_("Expected nrow(Ztl[[%d]]) to be %d, got %d"),
		  i + 1, Gp[i + 1] - Gp[i], ts1->nrow);
	ts2 = M_cholmod_vertcat(Zt, ts1, TRUE, &c);
	M_cholmod_free_sparse(&Zt, &c);
	Zt = ts2;
    }
				/* determine permutation */
    Perm = Alloca(q + p + 1, int);
    R_CheckStack();
    for (j = 0; j <= (p + q); j++) Perm[j] = j; /* initialize to identity */
    Zt_perm(Gp, Zt, ST, dims, Perm);
				/* create [X;-y]' */
    Xy = AS_CHM_DN(Xp);
    R_CheckStack();
    ts1 = M_cholmod_dense_to_sparse(Xy, 1 /* values */, &c);
    ts2 = M_cholmod_transpose(ts1, 1 /* values */, &c);
    M_cholmod_free_sparse(&ts1, &c);
    ts1 = ts2;
    /* y positions in ZXyt are always stored, even if zero */
    ts2 = M_cholmod_allocate_sparse(ts1->nrow + 1, ts1->ncol,
				    ts1->nzmax + ts1->ncol, ts1->sorted,
				    ts1->packed, 0 /*stype */,
				    CHOLMOD_REAL, &c);
    p1 = (int*)(ts1->p);
    p2 = Memcpy((int*)(ts2->p), p1, n + 1);
    for (j = 0; j <= n; j++) p2[j] += j; /* 1 extra nz per col */
    i1 = (int*)(ts1->i); i2 = (int*)(ts2->i);
    src = (double*)(ts1->x); dest = (double*)(ts2->x);
    for (j = 0; j < n; j++) {
	for (i = p1[j]; i < p1[j + 1]; i++) {
	    *i2++ = i1[i];
	    *dest++ = src[i];
	}
	*i2++ = ts1->nrow;
	*dest++ = (hasY) ? -(*y++) : 0;
    }
    M_cholmod_free_sparse(&ts1, &c);
    ts1 = ts2;
    ts2 = M_cholmod_vertcat(Zt, ts1, TRUE, &c);
    M_cholmod_free_sparse(&Zt, &c);
    M_cholmod_free_sparse(&ts1, &c);
    SET_SLOT(val, lme4_ZXytSym,
	     M_chm_sparse_to_SEXP(ts2, 0, 0, 0, "", R_NilValue));
				/*  Create and store A's pattern */
    ts1 = M_cholmod_aat(ts2, (int*)NULL, (size_t)0, 1, &c);
    ts3 = M_cholmod_copy(ts1, +1/*upper triangle*/, +1/*values*/, &c);
    M_cholmod_free_sparse(&ts1, &c);
    make_cholmod_sparse_sorted(ts3);
    SET_SLOT(val, lme4_ASym,
	     M_chm_sparse_to_SEXP(ts3, 0, 0, 0, "", R_NilValue));
    M_cholmod_free_sparse(&ts3, &c);
				/* update A using weights and offset */
    A = AS_CHM_SP(GET_SLOT(val, lme4_ASym));
    R_CheckStack();
    internal_update_A(ts2, GET_SLOT(val, lme4_weightsSym),
		      GET_SLOT(val, lme4_offsetSym), A);
    M_cholmod_free_sparse(&ts2, &c);
    internal_lmer2_initial(ST, Gp, A);
    j = c.supernodal;		/* force supernodal for non-nested */
    if (!dims[isNest_POS]) c.supernodal = CHOLMOD_SUPERNODAL;
    i = c.nmethods;
    c.nmethods = 1;		/* force user-specified permutation */
    				/* Create L  with user-specified perm */
    L = M_cholmod_analyze_p(A, Perm, (int*)NULL, (size_t)0, &c);
    if (!L)
	error(_("cholmod_analyze_p returned with status %d"), c.status);
				/* restore previous settings */
    c.nmethods = i;
    c.supernodal = j;
				/* initialize and store L */
    internal_update_L(REAL(GET_SLOT(val, lme4_devianceSym)),
			  dims, Gp, ST, A, L);
    SET_SLOT(val, lme4_LSym, M_chm_factor_to_SEXP(L, 1));
}

/**
 * Check a numeric vector (of weights) to see if it is all ones
 *
 * @param x Pointer to a numeric vector
 *
 * @return x if it contains any values that are not 1, otherwise numeric(0)
 */
static SEXP all_ones(SEXP x)
{
    int i, n = LENGTH(x);
    double *xx = REAL(x);

    if (!isReal(x))
	error(_("argument x to all_ones must be numeric"));
    for (i = 0; i < n; i++) if (xx[i] != 1) return x;
    return allocVector(REALSXP, 0);
}


/**
 * Remove the names attribute from an object
 *
 * @param x 
 *
 * @return x without the names attribute
 */
static SEXP R_INLINE unname(SEXP x) {
    setAttrib(x, R_NamesSymbol, R_NilValue);
    return x;
}

/**
 * Create an lmer2 object from a list of grouping factors and a list of model
 * matrices.
 *
 * @param fr list produced by lmerFrames
 * @param FL factor list produced by lmerFactorList
 * @param glmp list produced by glm.fit
 * @param method character vector
 * @param mc matched call
 * @param mod logical scalar indicating if the model frame should be
 *        saved in the returned object.
 *
 * @return pointer to an lmer2 object
 */
SEXP lmer2_create(SEXP fr, SEXP FL, SEXP Ztl, SEXP glmp,
		  SEXP method, SEXP mc, SEXP mod)
{
    SEXP Xp = getListElement(fr, "X"),
	Ztorig = getListElement(FL, "Ztl"),
	family = getListElement(glmp, "family"),
	fl = getListElement(FL, "fl"),
	offset = getListElement(fr, "offset"),
	wts = getListElement(glmp, "prior.weights"),
	yp = getListElement(glmp, "y");
    int *xdims = INTEGER(getAttrib(Xp, R_DimSymbol)), *Gp, *dims,
	REML = !strcmp(CHAR(asChar(method)), "REML"),
	ftyp = flType(family), i, nf = LENGTH(fl), nobs = LENGTH(yp), p, q;
    SEXP ST, Xdnames, cnames, cnamesnames, flnms, glmFixed,
	val = PROTECT(NEW_OBJECT(MAKE_CLASS(ftyp < 0 ? "lmer2" : "glmer2")));
    char *DEVIANCE_NAMES[]={"ML","REML","ldZ","ldX","lr2", "bqd", "Sdr", ""};
    char *DIMS_NAMES[]={"nf","n","p","q", "s", "np","isREML","famType","isNested",""};
    double *awv /* adjusted working variate */, *irlsw /* IRLS wts */;
				/* check arguments */
    if (!isReal(yp) || nobs <= 0)
	error(_("y must be a non-null numeric vector"));
    if (!isNewList(fl) || nf <= 0)
	error(_("fl must be a non-null list"));
    if (!isMatrix(Xp) || (xdims[1] && !isReal(Xp)))
	error(_("X must be a numeric matrix"));
    if (*xdims != nobs)
	error(_("Dimension mismatch: length(y) = %d and nrow(X) = %d"),
	      nobs, xdims[0]);
    if (!isNewList(Ztl) || LENGTH(Ztl) != nf) 
	error(_("Length mismatch: fl has length %d and Ztl has length %d"),
	      nf, LENGTH(Ztl));
    if (!isNewList(Ztorig) || LENGTH(Ztorig) != nf)
	error(_("FL$Ztl must be a list of length %d"), nf); 
				/*  Copy arguments to slots*/
    SET_SLOT(val, install("frame"),
	     (asLogical(mod)) ? duplicate(getListElement(fr, "mf"))
	     : NullFrame());
    SET_SLOT(val, install("terms"), duplicate(getListElement(fr, "mt")));
    SET_SLOT(val, install("call"), duplicate(mc));
				/* allocate slots */
    flnms = getAttrib(fl, R_NamesSymbol);
    ST = ALLOC_SLOT(val, lme4_STSym, VECSXP, nf);
    setAttrib(ST, R_NamesSymbol, duplicate(flnms));
    /* set up the cnames - tedious but convenient to have */
    cnames = ALLOC_SLOT(val, lme4_cnamesSym, VECSXP, nf + 2);
    setAttrib(cnames, R_NamesSymbol, allocVector(STRSXP, nf + 1));
    cnamesnames = getAttrib(cnames, R_NamesSymbol);
    for (i = 0; i < nf; i++)
	SET_STRING_ELT(cnamesnames, i, STRING_ELT(flnms, i));
    SET_STRING_ELT(cnamesnames, nf, mkChar(".fixed"));
    Xdnames = getAttrib(Xp, R_DimNamesSymbol);
    if (!isNull(Xdnames) && isNewList(Xdnames) && LENGTH(Xdnames) == 2)
	SET_VECTOR_ELT(cnames, nf, duplicate(VECTOR_ELT(Xdnames, 1)));
				/* dims vector */
    SET_SLOT(val, lme4_dimsSym, internal_make_named(INTSXP, DIMS_NAMES));
    dims = INTEGER(GET_SLOT(val, lme4_dimsSym));
    dims[nf_POS] = nf;
    dims[n_POS] = nobs;
    dims[np_POS] = 0;
    dims[isREML_POS] = REML;
    dims[famType_POS] = ftyp;
    dims[isNest_POS] = TRUE;
				/* Create Gp and populate ST */
    Gp = INTEGER(ALLOC_SLOT(val, lme4_GpSym, INTSXP, nf + 3));
    Gp[0] = 0;
    for (i = 0; i < nf; i++) {
	SEXP fli = VECTOR_ELT(fl, i);
	SEXP Zti = VECTOR_ELT(Ztorig, i);
	SEXP Zdimnms = getAttrib(Zti, R_DimNamesSymbol);
	int *Zdims = INTEGER(getAttrib(Zti, R_DimSymbol));

	if (!isReal(Zti) || !isMatrix(Zti) || Zdims[1] != nobs || Zdims[0] <= 0)
	    error(_("FL$Ztl[[%d]] is not a numeric matrix with %n columns and > 0 rows"),
		  i + 1, nobs);
	SET_VECTOR_ELT(ST, i, allocMatrix(REALSXP, Zdims[0], Zdims[0]));
	dims[np_POS] += (Zdims[0] * (Zdims[0] + 1))/2;
	if (isNewList(Zdimnms) && LENGTH(Zdimnms) == 2)
	    SET_VECTOR_ELT(cnames, i, duplicate(VECTOR_ELT(Zdimnms, 0)));
	if (!isFactor(fli) || LENGTH(fli) != nobs)
	    error(_("fl[[%d] must be a factor of length %d"), i+1, nobs);
	Gp[i + 1] = Gp[i] + LENGTH(getAttrib(fli, R_LevelsSymbol)) * Zdims[0];
    }
    dims[p_POS] = p = (ftyp < 0) ? xdims[1] : 0;
    dims[q_POS] = q = Gp[nf];
    Gp[nf + 1] = Gp[nf] + p;	 /* fixed effects */
    Gp[nf + 2] = Gp[nf + 1] + 1; /* response */
    SET_SLOT(val, lme4_flistSym, duplicate(fl));

    AZERO(REAL(ALLOC_SLOT(val, lme4_ranefSym, REALSXP, q)), q);
    /* initialize fixed effects (not needed for lmm but harmless) */
    glmFixed = getListElement(glmp, "coefficients");
    if (LENGTH(glmFixed) != xdims[1])
	error(_("Dimension mismatch: length(coef) = %d != ncol(X) = %d"),
	      LENGTH(glmFixed), xdims[1]);
    SET_SLOT(val, lme4_fixefSym, unname(duplicate(glmFixed)));
    SET_SLOT(val, lme4_weightsSym,
	     (wts == R_NilValue) ? allocVector(REALSXP, 0) :
	     unname(duplicate(all_ones(wts))));
    SET_SLOT(val, lme4_offsetSym,
	     (offset == R_NilValue) ? allocVector(REALSXP, 0)
	     : unname(duplicate(offset)));
    SET_SLOT(val, lme4_devianceSym,
	     internal_make_named(REALSXP, DEVIANCE_NAMES));
    if (dims[famType_POS] < 0) { /* linear mixed model */
	ZXyt_create(Ztl, Xp, yp, val);
	UNPROTECT(1);
	return val;
    }
				/* generalized linear mixed model */
    dims[isREML_POS] = 0;	/* REML not defined for GLMMs */
    SET_SLOT(val, install("family"), duplicate(family));
    SET_SLOT(val, lme4_ySym, unname(duplicate(yp)));
    SET_SLOT(val, lme4_XSym, duplicate(Xp));
    SET_SLOT(val, install("moff"),  /* offset in model, if any */
	     GET_SLOT(val, lme4_offsetSym)); 
    SET_SLOT(val, install("pwts"), /* prior weights */
	     unname(duplicate(getListElement(glmp, "prior.weights"))));
    SET_SLOT(val, lme4_etaSym,
	     unname(duplicate(getListElement(glmp, "linear.predictors"))));
    SET_SLOT(val, lme4_muSym,
	     unname(duplicate(getListElement(glmp, "fitted.values"))));
    /* non-null values in adj. wrk. variate for creating A and L */
    awv = REAL(ALLOC_SLOT(val, lme4_offsetSym, REALSXP, nobs));
    /* non-null values in IRLS weights for creating A and L */
    irlsw = REAL(ALLOC_SLOT(val, lme4_weightsSym, REALSXP, nobs));
    for (i = 0; i < nobs; i++) {awv[i] = 0.5; irlsw[i] = 1.0;}

    /* Empty X and y slots when creating ZXyt (y's are stored as 0's) */
    ZXyt_create(Ztl, PROTECT(allocMatrix(REALSXP, nobs, 0)),
		PROTECT(allocVector(REALSXP, 0)), val);
    glmer_reweight(val);
    UNPROTECT(3);
    return val;
}

/**
 * Update the response in an lmer2 object.
 *
 * @param x an lmer2 object
 * @param y a numeric vector of length nobs to update y
 *
 * @return R_NilValue
 */
SEXP lmer2_update_y(SEXP x, SEXP yp)
{
    CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym)),
	ZXyt = AS_CHM_SP(GET_SLOT(x, lme4_ZXytSym));
    int *ip = (int*)(ZXyt->i), *pp = (int*)(ZXyt->p), j;
    int n = ZXyt->ncol, yrow = ZXyt->nrow - 1;
    double *xp = (double*)(ZXyt->x), *y = REAL(yp);
    R_CheckStack();

    if (!isReal(yp) || LENGTH(yp) != n)
	error(_("y must be a numeric vector of length %d"), n);
    for (j = 0; j < n; j++) {
	int ind = pp[j + 1] - 1;
	if (ip[ind] != yrow)
	    error(_("Missing y position in column %d of ZXyt"), j + 1);
	xp[ind] = -y[j];
    }
    internal_update_A(ZXyt, GET_SLOT(x, lme4_weightsSym),
		      GET_SLOT(x, lme4_offsetSym), A);
    return R_NilValue;
}

#define IRLS_MAXITER  60
#define IRLS_TOL      1e-9

/**
 * Evaluate the convergence criterion and copy eta to
 * etaold
 *
 * @param eta current value of the linear predictors
 * @param etaold previous values of the linear predictors
 *
 * @return convergence criterion
 */
static double
conv_crit(double etaold[], double eta[], int n) {
    double max_abs_eta = -1, max_abs_diff = -1;
    int i;

    for (i = 0; i < n; i++) {
	double abs_eta, abs_diff;

	abs_eta = fabs(eta[i]);
	if (abs_eta > max_abs_eta) max_abs_eta = abs_eta;
	abs_diff = fabs(eta[i] - etaold[i]);
	if (abs_diff > max_abs_diff) max_abs_diff = abs_diff;
	etaold[i] = eta[i];
    }
    return max_abs_diff / (0.1 + max_abs_eta);
}

/**
 * Iterate to determine the conditional modes of the random effects.
 *
 * @param x pointer to a glmer2 object
 *
 * @return An indicator of whether the iterations converged
 */
static int internal_bhat(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, n = dims[n_POS];
    double *dev = REAL(GET_SLOT(x, lme4_devianceSym)),
	*eta = REAL(GET_SLOT(x, lme4_etaSym)),
	*etaold = Calloc(n, double), crit = IRLS_TOL + 1;
#ifdef DEBUG_LMER2
    double *dev_res = Alloca(n, double);
#endif
    CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym));
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    R_CheckStack();

    glmer_eta(x);
    glmer_reweight(x);
    internal_update_L(dev, dims, Gp, ST, A, L);
    lmer2_update_effects(x);
    glmer_eta(x);
    Memcpy(etaold, eta, n);
    for (i = 0; i < IRLS_MAXITER && crit > IRLS_TOL; i++) {
	glmer_reweight(x);
	internal_update_L(dev, dims, Gp, ST, A, L);
	lmer2_update_effects(x);
	glmer_eta(x);
#ifdef DEBUG_LMER2
	glmer_linkinv(x);
	Rprintf("%3d: %20.15g %20.15g", i,
		dev[Sdr_POS] = glmer_dev_resids(x, dev_res), dev[bqd_POS]);
	Rprintf(" %20.15g\n", dev[Sdr_POS] + dev[bqd_POS]);
#endif
	crit = conv_crit(etaold, eta, n);
    }
    return (crit > IRLS_TOL) ? 0 : i;
}

/**
 * Determine the conditional modes of the random effects.
 *
 * @param x pointer to a glmer2 object
 *
 * @return Number of iterations to convergence (0 for non-convergence)
 */
SEXP glmer_bhat(SEXP x) {
    return ScalarInteger(internal_bhat(x));
}

/**
 * Factor (Vt WW (Vt)' + I) as L, solve
 * (LL')u = (Vt W z) where W is diagonal
 *
 * @param Vt - sparse representation of V'
 * @param w - diagonal of matrix W
 * @param z - working residual vector
 * @param L - factor to update
 * @param u - holds the solution
 *
 * @return squared length of u
 */
static double
Vt_update_L(cholmod_sparse *Vt, const double *w,
	     const double *z,
	     cholmod_factor *L, double *u)
{
    cholmod_sparse *wZ = M_cholmod_copy_sparse(Vt, &c);
    int *wi = (int*)(wZ->i), *wp = (int*)(wZ->p),
	i, j, m = wZ->nrow, n = wZ->ncol;
    cholmod_dense *td,
	*rhs = M_cholmod_allocate_dense(m, 1, m, CHOLMOD_REAL, &c);
    double *wx = (double*)(wZ->x), *rh = (double*)(rhs->x),
	one[] = {1, 0}, val;

    AZERO(rh, m);
    for (j = 0; j < n; j++) {
	for (i = wp[j]; i < wp[j + 1]; i++) {
	    wx[i] *= w[j];	/* weight jth column */
	    rh[wi[i]] += wx[i] * z[j]; /* inner product */
	}
    }
    if (!M_cholmod_factorize_p(wZ, one, (int*) NULL, (size_t) 0, L, &c)) { 
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    }
    M_cholmod_free_sparse(&wZ, &c);
				/* CHOLMOD_A also applies any permutation */
    td = M_cholmod_solve(CHOLMOD_A, L, rhs, &c);
    M_cholmod_free_dense(&rhs, &c);
    for (i = 0, val = 0; i < m; i++) {
	double tt = ((double*)(td->x))[i];
	u[i] = tt;
	val += tt * tt;
    }
    M_cholmod_free_dense(&td, &c);
    return val;
}

/**
 * Iterate to determine the conditional modes of the random effects.
 *
 * @param x pointer to a glmer2 object
 *
 * @return An indicator of whether the iterations converged
 */
static int internal_bhat2(SEXP x)
{
    CHM_SP Vt = AS_CHM_SP(GET_SLOT(x, install("Vt")));
    SEXP etap = GET_SLOT(x, lme4_etaSym),
	fixef = GET_SLOT(x, lme4_fixefSym),
	moff = GET_SLOT(x, install("moff")), /* this can revert to offset */
	ranef = GET_SLOT(x, lme4_ranefSym),
	wts = GET_SLOT(x, lme4_weightsSym);
    int i, ione = 1, n = Vt->ncol, p = LENGTH(fixef);
    double *Xbeta = Calloc(n, double),
	*dev = REAL(GET_SLOT(x, lme4_devianceSym)),
	*dmu_deta = Calloc(n, double),
	*eta = REAL(etap), 
	*mu = REAL(GET_SLOT(x, lme4_muSym)),
	*pwts = REAL(GET_SLOT(x, install("pwts"))), /* can revert to weights */
	*u = REAL(ranef),
	*var = Calloc(n, double), *w = REAL(wts),
	*y = REAL(GET_SLOT(x, lme4_ySym)), *z = Calloc(n, double),
	*etaold = Calloc(n, double), crit = IRLS_TOL + 1, one[] = {1,0};
    CHM_DN ceta = AS_CHM_DN(etap), cu = AS_CHM_DN(ranef);
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
#ifdef DEBUG_LMER2
    double *dev_res = Alloca(n, double);
#endif
    R_CheckStack();
				
    AZERO(Xbeta, n);		/* Evaluate offset = moff + X\beta */
    if (LENGTH(moff) == n) Memcpy(Xbeta, REAL(moff), n);
    F77_CALL(dgemv)("N", &n, &p, one, REAL(GET_SLOT(x, lme4_XSym)),
		    &n, REAL(fixef), &ione, one, Xbeta, &ione);

    Memcpy(eta, Xbeta, n);	/* initialize eta to Xbeta */
    Memcpy(etaold, eta, n);
    for (i = 0; i < IRLS_MAXITER && crit > IRLS_TOL; i++) {
	glmer_linkinv(x);	/* evaluate mu */
	glmer_dmu_deta(x, dmu_deta);
	glmer_var(x, var);
	Memcpy(w, pwts, n);
	for (i = 0; i < n; i++) {
	    w[i] = sqrt(pwts[i] * dmu_deta[i] * dmu_deta[i]/var[i]);
	    z[i] = w[i] * (eta[i] - Xbeta[i] + (y[i] - mu[i])/dmu_deta[i]);
	}
	dev[bqd_POS] = Vt_update_L(Vt, w, z, L, u);
	if (!M_cholmod_sdmult(Vt, 1 /* trans */, one, one, cu, ceta, &c))
	    error(_("cholmod_sdmult error returned"));
#ifdef DEBUG_LMER2
	glmer_linkinv(x);
	Rprintf("%3d: %20.15g %20.15g", i,
		dev[Sdr_POS] = glmer_dev_resids(x, dev_res), dev[bqd_POS]);
	Rprintf(" %20.15g\n", dev[Sdr_POS] + dev[bqd_POS]);
#endif
	crit = conv_crit(etaold, eta, n);
    }
    Free(Xbeta); Free(dmu_deta); Free(var); Free(etaold); Free(z);
/* FIXME: Should this be freed with cholmod_free_sparse?  I don't
   think it has been changed from the initial AS_CHM_SP allocation. */
/*     M_cholmod_free_sparse(&Vt, &c); */
    return (crit > IRLS_TOL) ? 0 : i;

}

/**
 * Determine the conditional modes of the random effects.
 *
 * @param x pointer to a glmer2 object
 *
 * @return Number of iterations to convergence (0 for non-convergence)
 */
SEXP glmer_bhat2(SEXP x) {
    return ScalarInteger(internal_bhat2(x));
}


/**
 * Evaluate starting estimates for the elements of ST
 *
 * @param ST pointers to the nf ST factorizations of the diagonal
 *     elements of Sigma 
 * @param Gp length nf+1 vector of group pointers for the rows of A
 * @param Zt transpose of Z matrix
 *
 */
static void
internal_nlmer_initial(SEXP ST, int *Gp, SEXP Zt)
{
    int	*Zdims = INTEGER(GET_SLOT(Zt, lme4_DimSym)),
	*zi = INTEGER(GET_SLOT(Zt, lme4_iSym)),
	i, nf = LENGTH(ST);
    int nnz = INTEGER(GET_SLOT(Zt, lme4_pSym))[Zdims[1]];
    double *rowsqr = Calloc(Zdims[0], double),
	*zx = REAL(GET_SLOT(Zt, lme4_xSym));
    
    AZERO(rowsqr, Zdims[0]);
    for (i = 0; i < nnz; i++) rowsqr[zi[i]] += zx[i] * zx[i];
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int bb = Gp[i], j, ncip1 = nci + 1;

	
	AZERO(st, nci * nci);
	for (j = bb; j < Gp[i + 1]; j++)
	    st[((j - bb) % nci) * ncip1] += rowsqr[j];
	for (j = 0; j < nci; j++)
	    st[j * ncip1] =
		sqrt(((double)((Gp[i+1]-Gp[i])/nci))/(0.375*st[j * ncip1]));
    }
    Free(rowsqr);
}

/**
 * Evaluate the nonlinear model, its gradient matrix and the residual
 * sum of squares
 *
 * @param x pointer to an nlmer object
 * @param uform logical indicator of whether to use V'u (as opposed to Z'b)
 *
 * @return the residual sum of squares
 */
static double
internal_nlmer_eval_model(SEXP x, int uform)
{
    SEXP gg, pnames = GET_SLOT(x, install("pnames")),
	rho = GET_SLOT(x, install("env")), vv;
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), *gdims;
    int i, n = dims[n_POS], s = dims[s_POS];
    double *y = REAL(GET_SLOT(x, lme4_ySym)), *mu,
	one[] = {1, 0}, sumsq, zero[] = {0, 0};
    CHM_SP MM = AS_CHM_SP(GET_SLOT(x, install("Xt")));
    CHM_DN eff = AS_CHM_DN(GET_SLOT(x, lme4_fixefSym)),
	Phi = M_cholmod_allocate_dense(n * s, 1, n * s, CHOLMOD_REAL, &c);
    R_CheckStack();
    
				/* Evaluate Phi */
    if (!(i = M_cholmod_sdmult(MM, 1 /*trans*/, one, zero, eff, Phi, &c)))
	error(_("cholmod_sdmult returned error code %d"), i);
    MM = AS_CHM_SP(GET_SLOT(x, uform ? install("Vt") : lme4_ZtSym));
    eff = AS_CHM_DN(GET_SLOT(x, uform ? install("uvec") : lme4_ranefSym));
    R_CheckStack();
    if (!(i = M_cholmod_sdmult(MM, 1 /*trans*/, one, one, eff, Phi, &c)))
	error(_("cholmod_sdmult returned error code %d"), i);

    /* distribute the parameters in the environment */
    for (i = 0; i < dims[s_POS]; i++) {
	vv = findVarInFrame(rho, install(CHAR(STRING_ELT(pnames, i))));
	if (!isReal(vv) || LENGTH(vv) != n)
	    error(_("Parameter %s in the environment must be a length %d numeric vector"),
		  CHAR(STRING_ELT(pnames, i)), n);
	Memcpy(REAL(vv), ((double*)(Phi->x)) + i * n, n);
    }
    vv = PROTECT(eval(GET_SLOT(x, install("model")), rho));
    if (!isReal(vv) || LENGTH(vv) != n)
	error(_("evaluated model is not a numeric vector of length %d"), n);
    gg = getAttrib(vv, install("gradient"));
    if (!isReal(gg) || !isMatrix(gg))
	error(_("gradient attribute of evaluated model must be a numeric matrix"));
    gdims = INTEGER(getAttrib(gg, R_DimSymbol));
    if (gdims[0] != n ||gdims[1] != s)
	error(_("gradient matrix must be of size %d by %d"), n, s);
    SET_SLOT(x, lme4_muSym, vv);
    mu = REAL(vv);
    for (i = 0, sumsq = 0; i < n; i++) {
	double res = y[i] - mu[i];
	sumsq += res * res;
    }
    UNPROTECT(1);
    return sumsq;
}

SEXP nlmer_eval_model(SEXP x, SEXP uform)
{
    return ScalarReal(internal_nlmer_eval_model(x, asLogical(uform)));
}

static int*
nz_col(int *nz, int j, int nf, const int *Gp,
       const int *nc, const int *zi, const int *zp)
{
    int i, ii, k, nextra, zrow;
    for (ii = zp[j]; ii < zp[j + 1]; ii++) {
	k = Gp_grp(zrow = zi[ii], nf, Gp);
	nextra = (zrow - Gp[k]) % nc[k];
	for (i = 0; i <= nextra; i++) nz[zrow - i] = 1;
    }
    return nz;
}

SEXP nlmer_create_Vt(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym),
	Zt = GET_SLOT(x, lme4_ZtSym),
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix")));
    SEXP Zdims = GET_SLOT(Zt, lme4_DimSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)), *adims = INTEGER(Zdims),
	*vi, *vp, *zi = INTEGER(GET_SLOT(Zt, lme4_iSym)),
	*zp = INTEGER(GET_SLOT(Zt, lme4_pSym)),
	i, j, nf = LENGTH(ST), nnz;
    int *nc = Alloca(nf, int), *nz = Alloca(adims[0], int);
    R_CheckStack();
    
    SET_SLOT(ans, lme4_DimSym, duplicate(Zdims));
    SET_SLOT(ans, lme4_DimNamesSym, allocVector(VECSXP, 2));
    for (i = 0; i < nf; i++) 	/* populate nc */
	nc[i] = *INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol));
				/* create and evaluate the p slot */
    vp = INTEGER(ALLOC_SLOT(ans, lme4_pSym, INTSXP, adims[1] + 1));
    vp[0] = 0;
    for (j = 0; j < adims[1]; j++) {
	AZERO(nz, adims[0]);
	nz_col(nz, j, nf, Gp, nc, zi, zp);
	for (i = 0, nnz = 0; i < adims[0]; i++)
	    if (nz[i]) nnz++;
	vp[j + 1] = vp[j] + nnz;
    }
    vi = INTEGER(ALLOC_SLOT(ans, lme4_iSym, INTSXP, vp[adims[1]]));
    AZERO(REAL(ALLOC_SLOT(ans, lme4_xSym, REALSXP, vp[adims[1]])), vp[adims[1]]);
				/* fill in the i slot */
    for (j = 0; j < adims[1]; j++) {
	int pos = vp[j];
	AZERO(nz, adims[0]);
	nz_col(nz, j, nf, Gp, nc, zi, zp);
	for (i = 0; i < adims[0]; i++) if (nz[i]) vi[pos++] = i;
    }

    UNPROTECT(1);
    return ans;
}

/**
 * Create an nlmer object
 *
 * @param env environment in which to evaluate the model
 * @param model nonlinear model expression as a call
 * @param frame model frame
 * @param pnames character vector of parameter names
 * @param call matched call to the R nlmer function
 * @param flist factor list
 * @param Xt transpose of fixed-effects model matrix
 * @param Zt transpose of random-effects model matrix
 * @param y response vector
 * @param weights weights vector (may have length 0)
 * @param cnames list of column names
 * @param Gp integer vector of group pointers
 * @param fixef numeric vector of fixed effects
 */
SEXP nlmer_create(SEXP env, SEXP model, SEXP frame, SEXP pnames,
		  SEXP call, SEXP flist, SEXP Xt, SEXP Zt, SEXP y,
		  SEXP weights, SEXP cnames, SEXP Gp, SEXP fixef)
{
    SEXP ST, ans = PROTECT(NEW_OBJECT(MAKE_CLASS("nlmer")));
    char *DEVIANCE_NAMES[]={"ML","REML","ldZ","ldX","lr2", "bqd", "Sdr", ""};
    char *DIMS_NAMES[]={"nf","n","p","q", "s", "np","isREML","famType","isNested",""};
    int *Gpp = INTEGER(Gp), *Zdims = INTEGER(GET_SLOT(Zt, lme4_DimSym)), *dims, i, iT;
    cholmod_sparse *cVt, *ts1, *ts2;
    cholmod_factor *L;
    double one[] = {1,0};

    SET_SLOT(ans, install("env"), env);
    SET_SLOT(ans, install("model"), model);
    SET_SLOT(ans, install("frame"), frame);
    SET_SLOT(ans, install("pnames"), pnames);
    SET_SLOT(ans, install("call"), call);
    SET_SLOT(ans, lme4_flistSym, flist);
    SET_SLOT(ans, install("Xt"), Xt);
    SET_SLOT(ans, lme4_ZtSym, Zt);
    SET_SLOT(ans, lme4_ySym, y);
    SET_SLOT(ans, lme4_weightsSym, weights);
    SET_SLOT(ans, lme4_cnamesSym, cnames);
    SET_SLOT(ans, lme4_GpSym, Gp);
    SET_SLOT(ans, lme4_fixefSym, fixef);
    SET_SLOT(ans, lme4_dimsSym,
	     internal_make_named(INTSXP, DIMS_NAMES));
    dims = INTEGER(GET_SLOT(ans, lme4_dimsSym));
    dims[nf_POS] = LENGTH(flist);
    dims[n_POS] = LENGTH(y);
    dims[np_POS] = 0;
    dims[p_POS] = LENGTH(fixef);
    dims[q_POS] = Zdims[0];
    dims[s_POS] = Zdims[1]/dims[n_POS];
    dims[isREML_POS] = FALSE;
    dims[famType_POS] = -1;
    dims[isNest_POS] = TRUE;
    SET_SLOT(ans, lme4_ranefSym, allocVector(REALSXP, Zdims[0]));
    AZERO(REAL(GET_SLOT(ans, lme4_ranefSym)), Zdims[0]);
    SET_SLOT(ans, install("uvec"), allocVector(REALSXP, Zdims[0]));
    AZERO(REAL(GET_SLOT(ans, install("uvec"))), Zdims[0]);
    internal_nlmer_eval_model(ans, 0); /* check the model evaluation */

    SET_SLOT(ans, lme4_devianceSym,
	     internal_make_named(REALSXP, DEVIANCE_NAMES));
    AZERO(REAL(GET_SLOT(ans, lme4_devianceSym)), 7);

    SET_SLOT(ans, lme4_STSym, allocVector(VECSXP, dims[nf_POS]));
    ST = GET_SLOT(ans, lme4_STSym);
    iT = TRUE;			/* is T the identity? */
    for (i = 0; i < dims[nf_POS]; i++) {
	int nci = (Gpp[i + 1] - Gpp[i])/
	    LENGTH(getAttrib(VECTOR_ELT(flist, i), R_LevelsSymbol));
	SET_VECTOR_ELT(ST, i, allocMatrix(REALSXP, nci, nci));
	dims[np_POS] += (nci*(nci + 1))/2;
	if (nci > 1) iT = FALSE;
    }
    internal_nlmer_initial(ST, Gpp, Zt); /* initialize ST */
    SET_SLOT(ans, install("Vt"), iT ? duplicate(Zt) : nlmer_create_Vt(ans));
    nlmer_update_Vt(ans);
    cVt = AS_CHM_SP(GET_SLOT(ans, install("Vt")));
    R_CheckStack();
    
    /* Create Mt beginning with s identity matrices concatenated horizontally */
    ts1 = M_cholmod_allocate_sparse((size_t) dims[n_POS],
				    (size_t) Zdims[1], (size_t) Zdims[1],
				    1/*sorted*/, 1/*packed*/,
				    0/*stype*/, CHOLMOD_REAL, &c);
    for (i = 0; i < Zdims[1]; i++) {
	((int*)(ts1->p))[i] = i;
	((int*)(ts1->i))[i] = i % dims[n_POS];
	((double*)(ts1->x))[i] = 1;
    }
    ((int*)(ts1->p))[Zdims[1]] = Zdims[1];
    ts2 = M_cholmod_transpose(ts1, TRUE/*values*/, &c);
    M_cholmod_free_sparse(&ts1, &c);
    ts1 = M_cholmod_ssmult(cVt, ts2, /* Create pattern for Mt */
			   0 /*stype*/, 1 /*values*/, 1 /*sorted*/, &c);
    M_cholmod_free_sparse(&ts2, &c);
    SET_SLOT(ans, install("Mt"),
	     M_chm_sparse_to_SEXP(ts1, 0, 0, 0, "", R_NilValue));
/* FIXME: Check for nesting in Mt here? */
    i = c.final_ll;
    c.final_ll = 1;
    L = M_cholmod_analyze(ts1, &c); /* Create pattern for L */
    if (!M_cholmod_factorize_p(ts1, one, (int*) NULL, 0, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    if (!M_cholmod_change_factor(CHOLMOD_REAL, 1 /* to_ll */,
				 L->is_super, 1 /* packed */,
				 1 /* monotonic */, L, &c))
	error(_("cholmod_change_factor failed"));
    c.final_ll = i;
    SET_SLOT(ans, lme4_LSym, M_chm_factor_to_SEXP(L, 0));
    M_cholmod_free_factor(&L, &c);

    M_cholmod_free_sparse(&ts1, &c);
    nlmer_update_Vt(ans);
    UNPROTECT(1);
    return ans;
}

/**
 * Update the transpose of M from the current gradient and Vt
 *
 * @param x pointer to an nlmer object
 *
 * @return R_NilValue
 */
SEXP nlmer_update_Mt(SEXP x)
{
    SEXP Mt = GET_SLOT(x, install("Mt")),
	Vt = GET_SLOT(x, install("Vt"));
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)),
	*mi = INTEGER(GET_SLOT(Mt, lme4_iSym)),
	*mp = INTEGER(GET_SLOT(Mt, lme4_pSym)),
	*vi = INTEGER(GET_SLOT(Vt, lme4_iSym)),
	*vp = INTEGER(GET_SLOT(Vt, lme4_pSym)), jv;
    double *grad = REAL(getAttrib(GET_SLOT(x, lme4_muSym),
				  install("gradient"))),
	*mx = REAL(GET_SLOT(Mt, lme4_xSym)),
	*vx = REAL(GET_SLOT(Vt, lme4_xSym));
    int n = dims[n_POS], s = dims[s_POS];

    AZERO(mx, mp[n]);
    for (jv = 0; jv < n * s; jv++) {
	int iv, jm = jv % n, im;
	for (iv = vp[jv]; iv < vp[jv + 1]; iv++) {
	    for (im = mp[jm]; im < mp[jm + 1]; im++)
		if (mi[im] == vi[iv]) break;
	    if (im == mp[jm + 1])
		error(_("Structure of Mt incompatible with Vt, jv = %d, iv = %d"),
		      jv, iv);
	    mx[im] += grad[jv] * vx[iv];
	}
    }
    return R_NilValue;
}

/**
 * Update the working residual as y - mu + M u
 *
 * @param x pointer to an nlmer object
 * @param wrkres array to hold the updated working residual
 *
 * @return wrkres pointer
 */
static
double *internal_nlmer_update_wrkres(SEXP x, double *wrkres)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), i;
    double *mu = REAL(GET_SLOT(x, lme4_muSym)),
	*y = REAL(GET_SLOT(x, lme4_ySym)), one[] = {1,0};
    CHM_SP Mt = AS_CHM_SP(GET_SLOT(x, install("Mt")));
    CHM_DN cwrk = N_AS_CHM_DN(wrkres, dims[n_POS]),
	cu = AS_CHM_DN(GET_SLOT(x, install("uvec")));
    R_CheckStack();

    Memcpy(wrkres, y, dims[n_POS]);
    if (!(i = M_cholmod_sdmult(Mt, 1 /* trans */, one, one, cu, cwrk, &c)))
	error(_("cholmod_sdmult returned error code %d"), i);
    for (i = 0; i < dims[n_POS]; i++) wrkres[i] -= mu[i];
    return wrkres;
}

/**
 * Externally callable function to return the working residual
 *
 * @param x pointer to an nlmer object
 *
 * @return working residual
 */
SEXP nlmer_update_wrkres(SEXP x)
{
    SEXP ans = PROTECT(allocVector(REALSXP, LENGTH(GET_SLOT(x, lme4_ySym))));

    internal_nlmer_update_wrkres(x, REAL(ans));
    UNPROTECT(1);
    return ans;
}

/**
 * Iterate to determine the conditional modes of the random effects.
 *
 * @param x pointer to an nlmer object
 *
 * @return An indicator of whether the iterations converged
 */
static int internal_nbhat(SEXP x)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), i, j;
    int n = dims[n_POS], q = dims[q_POS], zq[] = {0, dims[q_POS]};
    double *dev = REAL(GET_SLOT(x, lme4_devianceSym)),
	*u = REAL(GET_SLOT(x, install("uvec"))),
	*uold = Alloca(q, double), *z = Alloca(n, double),
	crit = IRLS_TOL + 1, dn = (double)n, one[] = {1,0}, zero[] = {0,0};
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    CHM_DN cz = N_AS_CHM_DN(z, n), cu,
	cMtz = M_cholmod_allocate_dense(q, 1, q, CHOLMOD_REAL, &c);
    CHM_SP Mt = AS_CHM_SP(GET_SLOT(x, install("Mt")));
    R_CheckStack();

    nlmer_update_Vt(x);
    Memcpy(uold, u, q);
    for (i = 0; i < IRLS_MAXITER && crit > IRLS_TOL; i++) {
	dev[Sdr_POS] = internal_nlmer_eval_model(x, 1);
	for (j = 0, dev[bqd_POS] = 0; j < q; j++) dev[bqd_POS] += u[j] * u[j];
#ifdef DEBUG_NLMER
	Rprintf("%3d: %20.15g %20.15g %20.15g\n", i, dev[Sdr_POS], dev[bqd_POS],
		dev[Sdr_POS] + dev[bqd_POS]);
#endif
	nlmer_update_Mt(x);
	j = c.final_ll;
	c.final_ll = L->is_ll;
	if (!M_cholmod_factorize_p(Mt, one, (int*) NULL, 0, L, &c))
	    error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
		  c.status, L->minor, L->n);
	c.final_ll = j;
	internal_nlmer_update_wrkres(x, z);
	if (!(j = M_cholmod_sdmult(Mt, 0 /* trans */, one, zero, cz, cMtz, &c)))
	    error(_("cholmod_sdmult returned error code %d"), j);
	if (!(cu = M_cholmod_solve(CHOLMOD_A, L, cMtz, &c)))
	    error(_("cholmod_solve (CHOLMOD_A) failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
	Memcpy(u, (double*)(cu->x), q);
	M_cholmod_free_dense(&cu, &c); 
/* FIXME: replace this by an orthogonality convergence criterion */
 	crit = conv_crit(uold, u, q);
    }
    dev[lr2_POS] = log(dev[Sdr_POS] + dev[bqd_POS]);
    chm_log_abs_det2(&(dev[ldZ_POS]), 1, zq, L);
    dev[ML_POS] = dev[REML_POS] =
	dev[ldZ_POS] + dn * (1. + dev[lr2_POS] + log(2. * PI / dn));

    M_cholmod_free_dense(&cMtz, &c);
    return (crit > IRLS_TOL) ? 0 : i;
}

SEXP nlmer_bhat(SEXP x) {
    return ScalarInteger(internal_nbhat(x));
}

/**
 * Set the parameters in an lmer2 object and evaluate the deviance of
 * an lmm or the Laplace approximation to the deviance of a glmm.
 *
 * @param x pointer to a glmer2 object
 * @param xv vector of parameter values
 * @param nfe number of fixed-effects parameters updated
 * @param mtype model type: 0 -> lmm, 1 -> nlmm, 2 -> glmm
 *
 * @return deviance
 */
static double
update_deviance(SEXP x, const double *xv, int nfe, int mtype)
{
    SEXP ST = GET_SLOT(x, lme4_STSym),
	fixefp = GET_SLOT(x, lme4_fixefSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int n = dims[n_POS];
    double *dev = REAL(GET_SLOT(x, lme4_devianceSym));

    internal_lmer2_setPars(xv, ST); /* common parameters */
    Memcpy(REAL(fixefp), xv + dims[np_POS], nfe);
    switch(mtype) {
    case 0: {			  /* linear mixed model */
	CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym));
	CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
	R_CheckStack();

	internal_update_L(dev, dims, INTEGER(GET_SLOT(x, lme4_GpSym)),
			  ST, A, L);
	break;
    }
    case 1: {			  /* nonlinear mixed model */
	internal_nbhat(x);
	break;
    }
    case 2: {
	double *dev_res = Alloca(n, double);
	R_CheckStack();

	Memcpy(REAL(fixefp), xv + dims[np_POS], nfe);
	if (!internal_bhat(x)) {
	    warning(_("IRLS iterations for b-hat did not converge"));
	}
	dev[Sdr_POS] = glmer_dev_resids(x, dev_res);
	dev[ML_POS] = dev[bqd_POS] + dev[Sdr_POS] + dev[ldZ_POS];
	break;
    }
    default:
	error(_("Unknown form of model for update_deviance"));
    }
    return dev[dims[isREML_POS] ? REML_POS : ML_POS];
}

/**
 * Set the parameters in an lmer2 object and evaluate the deviance of
 * an lmm or the Laplace approximation to the deviance of a glmm.
 *
 * @param x pointer to a glmer2 object
 * @param xv vector of parameter values
 *
 * @return deviance
 */
static int
internal_optimize(SEXP x, int verb, int nfe, int mtype)
{
    SEXP ST = GET_SLOT(x, lme4_STSym),
	fixefp = GET_SLOT(x, lme4_fixefSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), i, j,
	nf = length(ST), pos;
    int np = dims[np_POS];
    int nv = np + nfe;
    int liv = S_iv_length(OPT, nv), lv = S_v_length(OPT, nv);
    int *iv = Calloc(liv, int);
    double *b = Calloc(2 * nv, double), *d = Calloc(nv, double),
	*fixef = REAL(fixefp),
	*g = (double*)NULL, *h = (double*)NULL,
	*v = Calloc(lv, double),
	*xv = internal_lmer2_getPars(ST, Calloc(nv, double)),
	fx = R_PosInf;

    Memcpy(xv + np, fixef, nfe);
				/* initialize the state vectors v and iv */
    S_Rf_divset(OPT, iv, liv, lv, v);
    if (verb) iv[OUTLEV] = 1;
				/* set the bounds to plus/minus Infty  */
    for (i = 0; i < nv; i++) {
	b[2*i] = R_NegInf; b[2*i+1] = R_PosInf; d[i] = 1;
    }
				/* reset lower bounds on elements of S */
    for (i = 0, pos = 0; i < nf; i++) {
	int nc = *INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol));
	for (j = 0; j < nc; j++) b[pos + 2*j] = 0;
	pos += nc * (nc + 1);
    }
    S_nlminb_iterate(b, d, fx, g, h, iv, liv, lv, nv, v, xv);
    while (iv[0] == 1 || iv[0] == 2) {
	fx = update_deviance(x, xv, nfe, mtype); 
	S_nlminb_iterate(b, d, fx, g, h, iv, liv, lv, nv, v, xv);
    }
    i = iv[0];
    Free(iv); Free(v); Free(xv); Free(b); Free(d);
    return i;
}

SEXP lmer2_optimize(SEXP x, SEXP verb)
{
    return ScalarInteger(internal_optimize(x, asInteger(verb), 0, 0));
}

SEXP nlmer_optimize(SEXP x, SEXP verb)
{
    return ScalarInteger(internal_optimize(x, asInteger(verb),
					   LENGTH(GET_SLOT(x, lme4_fixefSym)),
					   1));
}

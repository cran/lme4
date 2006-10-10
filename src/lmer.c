#include "lmer.h"
#include <float.h>

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

/* Some of these utility functions should be moved to the R sources */

/**
 * Symmetrize a matrix by copying the strict upper triangle into the
 * lower triangle.
 *
 * @param a pointer to a matrix in Fortran storage mode
 * @param nc number of columns (and rows and leading dimension) in the matrix
 *
 * @return a, symmetrized
 */
static double*
internal_symmetrize(double *a, int nc)
{
    int i, j;

    for (i = 1; i < nc; i++)
	for (j = 0; j < i; j++)
	    a[i + j*nc] = a[j + i*nc];
    return a;
}

/**
 * Create a named vector of type TYP
 *
 * @param TYP a vector SEXP type (e.g. REALSXP)
 * @param names names of list elements with null string appended
 *
 * @return pointer to a named vector of type TYP
 */
static SEXP
internal_make_named(int TYP, char **names)
{
    SEXP ans, nms;
    int i, n;

    for (n = 0; strlen(names[n]) > 0; n++) {}
    ans = PROTECT(allocVector(TYP, n));
    nms = PROTECT(allocVector(STRSXP, n));
    for (i = 0; i < n; i++) SET_STRING_ELT(nms, i, mkChar(names[i]));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ans;
}

/**
 * Return the element of a given name from a named list
 *
 * @param list
 * @param nm name of desired element
 *
 * @return element of list with name nm
 */
static SEXP
internal_getElement(SEXP list, char *nm) {
    SEXP names = getAttrib(list, R_NamesSymbol);
    int i;

    for (i = 0; i < LENGTH(list); i++)
	if (!strcmp(CHAR(STRING_ELT(names, i)), nm))
	    return(VECTOR_ELT(list, i));
    return R_NilValue;
}

/**
 * Allocate a 3-dimensional array
 *
 * @param mode The R mode (e.g. INTSXP)
 * @param nrow number of rows
 * @param ncol number of columns
 * @param nface number of faces
 *
 * @return A 3-dimensional array of the indicated dimensions and mode
 */
static SEXP alloc3Darray(SEXPTYPE mode, int nrow, int ncol, int nface)
{
    SEXP s, t;
    int n;

    if (nrow < 0 || ncol < 0 || nface < 0)
	error(_("negative extents to 3D array"));
    if ((double)nrow * (double)ncol * (double)nface > INT_MAX)
	error(_("alloc3Darray: too many elements specified"));
    n = nrow * ncol * nface;
    PROTECT(s = allocVector(mode, n));
    PROTECT(t = allocVector(INTSXP, 3));
    INTEGER(t)[0] = nrow;
    INTEGER(t)[1] = ncol;
    INTEGER(t)[2] = nface;
    setAttrib(s, R_DimSymbol, t);
    UNPROTECT(2);
    return s;
}

/**
 * Simulate the Cholesky factor of a standardized Wishart variate with
 * dimension p and df degrees of freedom.
 *
 * @param df degrees of freedom
 * @param p dimension of the Wishart distribution
 * @param ans array of size p * p to hold the result
 *
 * @return ans
 */
static double*
std_rWishart_factor(double df, int p, double ans[])
{
    int i, j, pp1 = p + 1;

    if (df < (double) p || p <= 0)
	error("inconsistent degrees of freedom and dimension");
    for (j = 0; j < p; j++) {	/* jth column */
	ans[j * pp1] = sqrt(rchisq(df - (double) j));
	for (i = 0; i < j; i++) ans[i + j * p] = norm_rand();
    }
    return ans;
}


/* Internally used utilities */

#define flag_not_factored(x) *INTEGER(GET_SLOT(x, lme4_statusSym)) = 0

/**
 * Create the diagonal blocks of the variance-covariance matrix of the
 * random effects
 *
 * @param Linv - cholmod_sparse representation of L^{-1}
 * @param nf - number of grouping factors
 * @param Gp - group pointers
 * @param nc - number of columns per factor
 * @param bVar - list of 3-d arrays to be filled in
 * @param uplo - "U" or "L" for upper or lower triangle
 */
static void
Linv_to_bVar(cholmod_sparse *Linv, const int Gp[], const int nc[],
	     SEXP bVar, const char uplo[])
{
    int *Lii = (int*)(Linv->i), *Lip = (int*)(Linv->p), i, nf = LENGTH(bVar);
    double *Lix = (double*)(Linv->x), one[] = {1,0}, zero[] = {0,0};

    for (i = 0; i < nf; i++) {
	int *ind, j, nci = nc[i], maxnnz = 0;
	int ncisqr = nci * nci, nlev = (Gp[i + 1] - Gp[i])/nci;
	double *bVi = REAL(VECTOR_ELT(bVar, i)), *tmp;

	AZERO(bVi, nlev * ncisqr);
	for (j = 0; j < nlev; j++) {
	    int nzm = Lip[Gp[i] + (j + 1) * nci] - Lip[Gp[i] + j * nci];
	    if (nzm > maxnnz) maxnnz = nzm;
	}
	ind = Calloc(maxnnz, int);
	tmp = Calloc(maxnnz * nci, double);
	for (j = 0; j < nlev; j++) {
	    int jj, k, kk;
	    int *ap = Lip + Gp[i] + j * nci;
	    int nr = ap[1] - ap[0];

	    AZERO(tmp, maxnnz * nci);
	    Memcpy(ind, Lii + ap[0], nr);
	    Memcpy(tmp, Lix + ap[0], nr);
	    for (jj = 1; jj < nci; jj++) {
		for (k = ap[jj]; k < ap[jj + 1]; k++) {
		    int aik = Lii[k];
		    for (kk = 0; kk < nr; kk++) {
			if (aik == ind[kk]) {
			    tmp[kk + jj * maxnnz] = Lix[k];
			    aik = -1;
			    break;
			}
		    }
		    if (aik >= 0) {	/* did not find the row index */
			ind[nr] = aik;
			tmp[nr + jj * maxnnz] = Lix[k];
			nr++;
		    }
		}
	    }
	    F77_CALL(dsyrk)(uplo, "T", &nci, &nr, one, tmp, &maxnnz,
			    zero, bVi + j * ncisqr, &nci);
	}
	Free(ind); Free(tmp);
    }
}

/**
 * Re-evaluate the factorization of the components of Omega and the
 * factorization of the crossproduct matrix - usually after an update.
 *
 * @param x pointer to an mer object
 */
/* NOTE: The factorization of the Omega components should be made
 * part of mer_factor.  Add an argument to mer_factor called force. */
static void
internal_mer_refactor(SEXP x)
{
    SEXP Omega = GET_SLOT(x, lme4_OmegaSym);
    int *nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	i, info, nf = LENGTH(GET_SLOT(x, lme4_flistSym));
    for (i = 0; i < nf; i++) {
	SEXP omgi = VECTOR_ELT(Omega, i);
	int nci = nc[i];
	double *choli;

/* NOTE: This operation depends on dpoMatrix_chol returning the
   Cholesky component of the factor slot, not a duplicate of it */
	choli = REAL(GET_SLOT(M_dpoMatrix_chol(omgi), lme4_xSym));
	Memcpy(choli, REAL(GET_SLOT(omgi, lme4_xSym)), nci * nci);
	F77_CALL(dpotrf)("U", &nci, choli, &nci, &info);
	/* Yes, you really do need to do that decomposition.
	   The contents of choli before the decomposition can be
	   from the previous value of Omegai. */
	if (info)
	    error(_("Omega[[%d]] is not positive definite"), i + 1);
    }
    flag_not_factored(x);
    mer_factor(x);
}

static void
internal_mer_RZXinv(SEXP x)
{
    SEXP RZXP = GET_SLOT(x, lme4_RZXSym);
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    cholmod_dense *RZX = M_as_cholmod_dense(RZXP), *tmp1;
    int *dims = INTEGER(GET_SLOT(RZXP, lme4_DimSym)),
	*Perm = (int *)(L->Perm);
    int i, j, p = dims[1], q = dims[0];
    int *iperm = Calloc(q, int);
    double *RZXinv = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXinvSym),
				   lme4_xSym)), m1[2] = {-1, 0};

				/* create the inverse permutation */
    for (j = 0; j < q; j++) iperm[Perm[j]] = j;
				/* solve system in L' */
    tmp1 = M_cholmod_solve(CHOLMOD_Lt, L, RZX, &c);
    /* copy columns of tmp1 to RZXinv applying the inverse permutation */
    for (j = 0; j < p; j++) {
	double *dest = RZXinv + j * q, *src = ((double*)(tmp1->x)) + j * q;
	for (i = 0; i < q; i++) dest[i] = src[iperm[i]];
    }
    M_cholmod_free_dense(&tmp1, &c);
    F77_CALL(dtrsm)("R", "U", "N", "N", &q, &p, m1,
		    REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym), lme4_xSym)),
		    &p, RZXinv, &q);
    Free(iperm); Free(RZX); Free(L);
}

static void
internal_mer_bVar(SEXP x)
{
    int q = LENGTH(GET_SLOT(x, lme4_rZySym));
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    cholmod_sparse *eye = M_cholmod_speye(q, q, CHOLMOD_REAL, &c), *Linv;
    int *Perm = (int *)(L->Perm), *iperm = Calloc(q, int), i;
				/* create the inverse permutation */
    for (i = 0; i < q; i++) iperm[Perm[i]] = i;
				/* apply iperm to the identity matrix */
    for (i = 0; i < q; i++) ((int*)(eye->i))[i] = iperm[i];
				/* Create Linv */
    Linv = M_cholmod_spsolve(CHOLMOD_L, L, eye, &c);
    M_cholmod_free_sparse(&eye, &c);
    Linv_to_bVar(Linv, INTEGER(GET_SLOT(x, lme4_GpSym)),
		 INTEGER(GET_SLOT(x, lme4_ncSym)),
		 GET_SLOT(x, lme4_bVarSym), "U");
    M_cholmod_free_sparse(&Linv, &c);
    Free(L); Free(iperm);
}

/**
 * Evaluate the quadratic form in b defined by Omega
 *
 * @param b vector of random effects
 * @param Omega - list of dpoMatrix objects defining the pattern for Omega
 * @param nf - number of grouping factors
 * @param Gp - group pointers
 * @param nc - number of columns per factor
 *
 * @return value of the quadratic form in b
 */
static double
b_quadratic(const double b[], SEXP Omega, const int Gp[], const int nc[])
{
    int i, ione = 1, nf = LENGTH(Omega);
    double ans = 0., one[] = {1.,0.};

    for (i = 0; i < nf; i++) {
	int nci = nc[i], ntot = Gp[i + 1] - Gp[i];
	int nlev = ntot/nci;
	double *bcp = Memcpy(Calloc(ntot, double), b + Gp[i], ntot),
	    *omgf = REAL(GET_SLOT(M_dpoMatrix_chol(VECTOR_ELT(Omega, i)),
				  lme4_xSym));

	F77_CALL(dtrmm)("L", "U", "N", "N", &nci, &nlev, one, omgf,
			&nci, bcp, &nci);
	ans += F77_CALL(ddot)(&ntot, bcp, &ione, bcp, &ione);
	Free(bcp);
    }
    return ans;
}

/**
 * Calculate the deviance for a linear mixed model at arbitrary
 * parameter values.  REML is not used - I don't think it applies.
 *
 * @param x a mixed-effects model representation
 * @param sigma standard deviation of per-observation noise
 * @param beta fixed-effects parameter vector
 *
 * @return deviance
 */
static
double lmm_deviance(SEXP x, double sigma, const double beta[])
{
    SEXP rXy = GET_SLOT(x, lme4_rXySym);
    int i, ione = 1, p = LENGTH(rXy);
    double *dcmp = REAL(GET_SLOT(x, lme4_devCompSym)),
	*betacp = Memcpy(Calloc(p, double), beta, p),
	*rXyp = REAL(rXy),
	sprss; /* scaled penalized rss */

    mer_factor(x);
    F77_CALL(dtrmv)("U", "N", "N", &p,
		    REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym), lme4_xSym)),
		    &p, betacp, &ione);
    sprss = exp(dcmp[3])/(sigma * sigma);
    for (i = 0; i < p; i++) {
	double ri = (rXyp[i] - betacp[i])/sigma;
	sprss += ri * ri;
    }
    Free(betacp);

    return dcmp[4] - dcmp[5] + dcmp[0] * log(2.*PI*sigma*sigma) + sprss;
}

/**
 * Calculate the length of the parameter vector (historically called "coef"
 * even though these are not coefficients).
 *
 * @param nf number of factors
 * @param nc number of columns in the model matrices for each factor
 *
 * @return total length of the coefficient vector
 */
static R_INLINE
int coef_length(int nf, const int nc[])
{
    int i, ans = 0;
    for (i = 0; i < nf; i++) ans += (nc[i] * (nc[i] + 1))/2;
    return ans;
}

/**
 * Find a variable of a given name in a given environment and check
 * that its length and mode are correct.
 *
 * @param rho Environment in which to find the variable
 * @param nm Name of the variable to find
 * @param mode Desired mode
 * @param len Desired length (ignored if <= 0)
 *
 * @return
 */
static
SEXP find_and_check(SEXP rho, SEXP nm, SEXPTYPE mode, int len)
{
    SEXP ans;
    if (R_NilValue == PROTECT(ans = findVarInFrame(rho, nm)))
	error(_("environment `rho' must contain an object `%s'"),
	      CHAR(PRINTNAME(nm)));
    if (TYPEOF(ans) != mode)
	error(_("object `%s' of incorrect type"),
	      CHAR(PRINTNAME(nm)));
    if (len > 0 && LENGTH(ans) != len)
	error(_("object `%s' must be of length %d"),
	      CHAR(PRINTNAME(nm)), len);
    UNPROTECT(1);
    return ans;
}

/**
 * Evaluate an expression in an environment, check that the length and
 * mode are as expected and store the result.
 *
 * @param fcn expression to evaluate
 * @param rho environment in which to evaluate it
 * @param vv position to store the result
 *
 * @return vv with new contents
 */
static
SEXP eval_check_store(SEXP fcn, SEXP rho, SEXP vv)
{
    SEXP v = PROTECT(eval(fcn, rho));
    if (TYPEOF(v) != TYPEOF(vv) || LENGTH(v) != LENGTH(vv))
	error(_("fcn produced mode %d, length %d - wanted mode %d, length %d"),
	      TYPEOF(v), LENGTH(v), TYPEOF(vv), LENGTH(vv));
    switch (TYPEOF(v)) {
    case LGLSXP:
	Memcpy(LOGICAL(vv), LOGICAL(v), LENGTH(vv));
	break;
    case INTSXP:
	Memcpy(INTEGER(vv), INTEGER(v), LENGTH(vv));
	break;
    case REALSXP:
	Memcpy(REAL(vv), REAL(v), LENGTH(vv));
	break;
    default:
	error(_("invalid type for eval_check_store"));
    }
    UNPROTECT(1);
    return vv;
}

/**
 * Evaluate an expression in an environment, check that the length and
 * mode are as expected and return the result.
 *
 * @param fcn expression to evaluate
 * @param rho environment in which to evaluate it
 * @param mode desired mode
 * @param len desired length
 *
 * @return evaluated expression
 */
static SEXP
eval_check(SEXP fcn, SEXP rho, SEXPTYPE mode, int len) {
    SEXP v = PROTECT(eval(fcn, rho));
    if (TYPEOF(v) != mode || LENGTH(v) != len)
	error(_("fcn produced mode %d, length %d - wanted mode %d, length %d"),
	      TYPEOF(v), LENGTH(v), mode, len);
    UNPROTECT(1);
    return v;
}

typedef struct glmer_struct
{
    SEXP cv;         /* control values */
    SEXP mer;	     /* mixed-effects representation */
    SEXP rho;        /* environment in which to evaluate the calls */
    SEXP eta;        /* linear predictor */
    SEXP mu;         /* mean vector */
    SEXP linkinv;    /* expression for inverse link evaluation */
    SEXP mu_eta;     /* expression for dmu/deta evaluation */
    SEXP LMEopt;     /* expression for LME optimization */
    SEXP dev_resids; /* expression for deviance residuals */
    SEXP var;        /* expression for variance evaluation */
    double *offset;  /* offset for GLM */
    double *wts;     /* prior weights for GLM */
    double *y;       /* copy of response vector */
    double *etaold;  /* previous value of eta for evaluating  */
    int n;	     /* length of the response vector */
    int p;	     /* length of fixed effects vector */
    int nf;	     /* number of grouping factors */
    int npar;        /* total length of the parameter */
    int niterEM;     /* default number of ECME iterations */
    int EMverbose;   /* logical indicator */
    int maxiter;     /* maximum number of IRLS iterations */
    double tol;      /* convergence tolerance for IRLS iterations */
} glmer_struct, *GlmerStruct;


/**
 * Calculate fitted values for the current fixed and random effects.
 *
 * @param x pointer to an mer object
 * @param initial initial values (i.e. an offset) or (double *) NULL
 * @param val array to hold the fitted values
 *
 * @return pointer to a numeric array of fitted values
 */
static double *
internal_mer_fitted(SEXP x, const double initial[], double val[])
{
    SEXP fixef = GET_SLOT(x, lme4_fixefSym),
	ranef = GET_SLOT(x, lme4_ranefSym);
    int ione = 1, n = LENGTH(GET_SLOT(x, lme4_ySym)), p = LENGTH(fixef);
    double *X = REAL(GET_SLOT(x, lme4_XSym)), one[] = {1,0};
    cholmod_sparse *Zt = M_as_cholmod_sparse(GET_SLOT(x, lme4_ZtSym));
    cholmod_dense *chv = M_numeric_as_chm_dense(val, n),
	*chb = M_as_cholmod_dense(ranef);

    mer_secondary(x);
    if (initial) Memcpy(val, initial, n);
    else AZERO(val, n);
    F77_CALL(dgemv)("N", &n, &p, one, X, &n, REAL(fixef), &ione, one, val, &ione);
    if (!M_cholmod_sdmult(Zt, 1, one, one, chb, chv, &c))
	error(_("Error return from cholmod_sdmult"));
    Free(chv); Free(chb); Free(Zt);
    return val;
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
static int
internal_mer_isNested(int nf, const int nc[], const int Gp[], const int p[])
{
    int **cnz = Calloc(nf, int*), ans = 1, i, j, k, nct;

    for (i = 0, nct = 0; i < nf; i++) { /* total number of columns */
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
    Free(cnz);
    return ans;
}

/**
 * Extract the coefficients
 *
 * @param x pointer to an mer object
 * @param ptyp parameter type to extract
 * @param ans vector to hold the extracted values
 *
 * @return ans
 */
static double *
internal_mer_coef(SEXP x, int ptyp, double ans[])
{
    SEXP Omega = GET_SLOT(x, lme4_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	i, nf = length(Omega), vind;

    vind = 0;			/* index in ans */
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1;
	if (nci == 1) {
	    double dd = REAL(GET_SLOT(VECTOR_ELT(Omega, i), lme4_xSym))[0];
	    ans[vind++] = ptyp ? ((ptyp == 1) ? log(dd) : 1./dd) : dd;
	} else {
	    if (ptyp) {	/* L log(D) L' factor of Omega[,,i] */
		int j, k, ncisq = nci * nci;
	        double *tmp = Memcpy(Calloc(ncisq, double),
				     REAL(GET_SLOT(M_dpoMatrix_chol(VECTOR_ELT(Omega, i)),
						   lme4_xSym)), ncisq);
		for (j = 0; j < nci; j++) {
		    double diagj = tmp[j * ncip1];
		    ans[vind++] = (ptyp == 1) ? (2. * log(diagj)) :
			1./(diagj * diagj);
		    for (k = j + 1; k < nci; k++) {
			tmp[j + k * nci] /= diagj;
		    }
		}
		for (j = 0; j < nci; j++) {
		    for (k = j + 1; k < nci; k++) {
			ans[vind++] = tmp[j + k * nci];
		    }
		}
		Free(tmp);
	    } else {		/* upper triangle of Omega[,,i] */
		int j, k, odind = vind + nci;
		double *omgi = REAL(GET_SLOT(VECTOR_ELT(Omega, i), lme4_xSym));

		for (j = 0; j < nci; j++) {
		    ans[vind++] = omgi[j * ncip1];
		    for (k = j + 1; k < nci; k++) {
			ans[odind++] = omgi[k*nci + j];
		    }
		}
		vind = odind;
	    }
	}
    }
    return ans;
}

/**
 * Set the coefficient vector and perform a factorization
 *
 * @param x pointer to an mer object
 * @param cc vector of coefficients to assign
 * @param ptyp indicator of the parametrization being used
 */
static
void internal_mer_coefGets(SEXP x, const double cc[], int ptyp)
{
    SEXP Omega = GET_SLOT(x, lme4_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	cind, i, nf = length(Omega);

    cind = 0;
    for (i = 0; i < nf; i++) {
	SEXP Omegai = VECTOR_ELT(Omega, i);
	int j, k, nci = nc[i], ncisq = nc[i] * nc[i];
	double *omgi = REAL(GET_SLOT(Omegai, lme4_xSym));

	if (nci == 1) {
	    double dd = cc[cind++];
	    *omgi = ptyp ? ((ptyp == 1) ? exp(dd) : 1./dd) : dd;
	} else {
	    int odind = cind + nci, /* off-diagonal index */
		ncip1 = nci + 1;

	    if (ptyp) {
		/* FIXME: Write this as an LDL decomposition */
		double *tmp = Calloc(ncisq, double),
		    diagj, one = 1., zero = 0.;

		AZERO(omgi, ncisq);
		for (j = 0; j < nci; j++) {
		    double dd = cc[cind++];
		    tmp[j * ncip1] = diagj =
			(ptyp == 1) ? exp(dd/2.) : sqrt(1./dd);
		    for (k = j + 1; k < nci; k++) {
			tmp[k*nci + j] = cc[odind++] * diagj;
		    }
		}
		F77_CALL(dsyrk)("U", "T", &nci, &nci, &one,
				tmp, &nci, &zero, omgi, &nci);
		Free(tmp);
	    } else {
		for (j = 0; j < nci; j++) {
		    omgi[j * ncip1] = cc[cind++];
		    for (k = j + 1; k < nci; k++) {
			omgi[k*nci + j] = cc[odind++];
		    }
		}
	    }
	    cind = odind;
	}
    }
    internal_mer_refactor(x);
}

/**
 * Evaluate current estimate of sigma from an mer object
 *
 * @param x pointer to an mer object
 * @param REML indicator of whether to use REML.
 *           < 0  -> determine REML or ML from x@status
 *           == 0 -> use ML unconditionally
 *           > 0  -> use REML unconditionally
 *
 * @return
 */
static double
internal_mer_sigma(SEXP x, int REML)
{
    double *dcmp = REAL(GET_SLOT(x, lme4_devCompSym));

    if (REML < 0)		/* get REML from x */
	REML = INTEGER(GET_SLOT(x, lme4_statusSym))[1];
    mer_factor(x);
    return exp(dcmp[3]/2)/sqrt(dcmp[0] - (REML ? dcmp[1] : 0));
}

/**
 * Update the derived quantities (ZtZ, ZtX, XtX, Zty, Xty
 * and dcmp[2] = y'y when Z, X, y, wts or wrkres has been changed.
 *
 * @param x pointer to an mer object
 * @param perm permutation from the Cholesky factor
 */

static void
internal_mer_update_ZXy(SEXP x, int *perm)
{
    SEXP Xp = GET_SLOT(x, lme4_XSym), ZtZ = GET_SLOT(x, lme4_ZtZSym),
	Ztyp = GET_SLOT(x, lme4_ZtySym);
    SEXP ZtZx = GET_SLOT(ZtZ, lme4_xSym),
	ZtZp = GET_SLOT(ZtZ, lme4_pSym), ZtZi = GET_SLOT(ZtZ, lme4_iSym);
    int *dims = INTEGER(getAttrib(Xp, R_DimSymbol)), i, ione = 1, j;
    int n = dims[0], nnz, p = dims[1], q = LENGTH(Ztyp);
    cholmod_sparse *ts1, *ts2,
	*Zt = M_as_cholmod_sparse(GET_SLOT(x, lme4_ZtSym));
    cholmod_sparse *Ztcp = M_cholmod_copy_sparse(Zt, &c);
    int *Zp = (int*)Ztcp->p;
    double *XtX = REAL(GET_SLOT(GET_SLOT(x, lme4_XtXSym), lme4_xSym)),
	*Xty = REAL(GET_SLOT(x, lme4_XtySym)),
	*ZtX = REAL(GET_SLOT(GET_SLOT(x, lme4_ZtXSym), lme4_xSym)),
	*Zty = REAL(Ztyp),
	*wts = REAL(GET_SLOT(x, lme4_wtsSym)),
	one[] = {1, 0}, zero[] = {0,0};
    cholmod_dense *td1, *Xd = M_as_cholmod_dense(Xp),
	*wkrd = M_as_cholmod_dense(GET_SLOT(x, lme4_wrkresSym));
    cholmod_dense *Xcp = M_cholmod_copy_dense(Xd, &c),
	*wkrcp = M_cholmod_copy_dense(wkrd, &c);
    double *X = (double*)(Xcp->x), *Ztx = (double*)(Ztcp->x),
	*wtres = (double*)(wkrcp->x);
				/* Apply weights */
    for (i = 0; i < n; i++)
	wtres[i] = ((double*)(wkrd->x))[i] * wts[i];
    for (j = 0; j < p; j++)
	for (i = 0; i < n; i++)
	    X[i + j * n] = ((double*)(Xd->x))[i + j * n] * wts[i];
    for (j = 0; j < n; j++)
	for (i = Zp[j]; i < Zp[j + 1]; i++)
	    Ztx[i] = ((double*)(Zt->x))[i] * wts[j];
    Free(Zt); Free(Xd); Free(wkrd);

				/* y'y */
    REAL(GET_SLOT(x, lme4_devCompSym))[2] =
	F77_CALL(ddot)(&n, wtres, &ione, wtres, &ione);
				/* ZtZ */
    ts1 = M_cholmod_aat(Ztcp, (int *) NULL, (size_t) 0, 1/* mode */, &c);
    /* cholmod_aat returns stype == 0; copy to set stype == 1 */
    ts2 = M_cholmod_copy(ts1, 1/* stype */, 1/* mode */, &c);
    nnz = M_cholmod_nnz(ts2, &c);
    if (((int)(ts2->ncol) + 1) != LENGTH(ZtZp))
	error(_("Order of Z'Z has changed - was %d, now %d"),
	      LENGTH(ZtZp) - 1, (int)(ts2->ncol));
    Memcpy(INTEGER(ZtZp), (int*)(ts2->p), LENGTH(ZtZp));
    if (nnz != LENGTH(ZtZx))
	error(_("Number of nonzeros in Z'Z has changed - was %d, now %d"),
	      LENGTH(ZtZx), nnz);
    Memcpy(INTEGER(ZtZi), (int*)(ts2->i), nnz);
    Memcpy(REAL(ZtZx), (double*)(ts2->x), nnz);
    M_cholmod_free_sparse(&ts1, &c); M_cholmod_free_sparse(&ts2, &c);
				/* PZ'X into ZtX */
    td1 = M_cholmod_allocate_dense(q, p, q, CHOLMOD_REAL, &c);
    if (!M_cholmod_sdmult(Ztcp, 0, one, zero, Xcp, td1, &c))
	error(_("Error return from cholmod_sdmult"));
    for (j = 0; j < p; j++) { 	/* apply the permutation to each column */
	double *dcol = ZtX + j * q,
	    *scol = (double*)(td1->x) + j * q;
	for (i = 0; i < q; i++) dcol[i] = scol[perm[i]];
    }
    M_cholmod_free_dense(&td1, &c);
				/* PZ'y into Zty */
    td1 = M_cholmod_allocate_dense(q, 1, q, CHOLMOD_REAL, &c);
    if (!M_cholmod_sdmult(Ztcp, 0, one, zero, wkrcp, td1, &c))
	error(_("Error return from cholmod_sdmult"));
    for (i = 0; i < q; i++) Zty[i] = ((double *)(td1->x))[perm[i]];
    M_cholmod_free_dense(&td1, &c); M_cholmod_free_sparse(&Ztcp, &c);
				/* XtX and Xty */
    AZERO(XtX, p * p);
    F77_CALL(dsyrk)("U", "T", &p, &n, one, X, &n, zero, XtX, &p);
    F77_CALL(dgemv)("T", &n, &p, one, X, &n, wtres, &ione, zero, Xty, &ione);
    M_cholmod_free_dense(&Xcp, &c); M_cholmod_free_dense(&wkrcp, &c);
    flag_not_factored(x);
}

static double chm_log_abs_det(cholmod_factor *F)
{
    double ans = 0;

    if (F->is_super) {
	int i;
	for (i = 0; i < F->nsuper; i++) {
	    int j, nrp1 = 1 + ((int *)(F->pi))[i + 1] - ((int *)(F->pi))[i],
		nc = ((int *)(F->super))[i + 1] - ((int *)(F->super))[i];
	    double *x = (double *)(F->x) + ((int *)(F->px))[i];

	    for (j = 0; j < nc; j++) ans += log(fabs(x[j * nrp1]));
	}
    } else
	error(_("code for simplicial factorization not yet written"));
    return ans;
}

static double
Omega_log_det(SEXP Omega, int *nc, int *Gp)
{
    double ans = 0;
    int i;

    for (i = 0; i < LENGTH(Omega); i++) {
	int j, nci = nc[i], ncip1 = nc[i] + 1, nlev = (Gp[i + 1] - Gp[i])/nc[i];
	double *omgi = REAL(GET_SLOT(M_dpoMatrix_chol(VECTOR_ELT(Omega, i)),
				     lme4_xSym));

	for (j = 0; j < nci; j++) ans += 2. * nlev * log(fabs(omgi[j * ncip1]));
    }
    return ans;
}


/**
 * Inflate Z'Z to Z'Z+Omega and factor. Form RZX and rZy and update
 * the status flags.
 *
 * @param x pointer to an mer object.
 */
static void
internal_mer_Zfactor(SEXP x, cholmod_factor *L)
{
    SEXP Omega = GET_SLOT(x, lme4_OmegaSym),
	rZyP = GET_SLOT(x, lme4_rZySym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	nf = LENGTH(Omega), q = LENGTH(rZyP),
	p = LENGTH(GET_SLOT(x, lme4_rXySym));
    cholmod_sparse *A, *Omg,
	*zz = M_as_cholmod_sparse(GET_SLOT(x, lme4_ZtZSym));
    cholmod_dense *ZtX = M_as_cholmod_dense(GET_SLOT(x, lme4_ZtXSym)),
	*Zty = M_as_cholmod_dense(GET_SLOT(x, lme4_ZtySym)),
	*rZy, *RZX;
    int *omp, *nnz = Calloc(nf + 1, int), i,
	*status = INTEGER(GET_SLOT(x, lme4_statusSym));
    double *dcmp = REAL(GET_SLOT(x, lme4_devCompSym)), one[] = {1, 0};


    dcmp[5] = Omega_log_det(Omega, nc, Gp); /* logDet(Omega) */
    for (nnz[0] = 0, i = 0; i < nf; i++)
	nnz[i + 1] = nnz[i] + ((Gp[i + 1] - Gp[i])*(nc[i] + 1))/2;
    Omg = M_cholmod_allocate_sparse(zz->nrow, zz->ncol, (size_t) nnz[nf],
				  TRUE, TRUE, 1, CHOLMOD_REAL, &c);
    omp = (int *) Omg->p;
    for (i = 0; i < nf; i++) {
	int bb = Gp[i], j, jj, k, nci = nc[i];
	int nlev = (Gp[i + 1] - bb)/nci;
	double *Omgi = REAL(GET_SLOT(VECTOR_ELT(Omega, i), lme4_xSym));

	for (j = 0; j < nlev; j++) { /* column of result */
	    int col0 = bb + j * nci; /* absolute column number */

	    for (jj = 0; jj < nci; jj++) { /* column of Omega_i */
		int coljj = col0 + jj;

		omp[coljj + 1] = omp[coljj] + jj + 1;
		for (k = 0; k <= jj; k++) { /* row of Omega_i */
		    int ind = omp[coljj];
		    ((int *)Omg->i)[ind + k] = col0 + k;
		    ((double *)Omg->x)[ind + k] = Omgi[jj * nci + k];
		}
	    }
	}
    }
    Free(nnz);
    A = M_cholmod_add(zz, Omg, one, one, TRUE, TRUE, &c);
    Free(zz); M_cholmod_free_sparse(&Omg, &c);
    if (!M_cholmod_factorize(A, L, &c))
	error(_("rank_deficient Z'Z+Omega"));
    M_cholmod_free_sparse(&A, &c);
    dcmp[4] = 2 * chm_log_abs_det(L); /* 2 * logDet(L) */

				/* calculate and store RZX and rZy */
    RZX = M_cholmod_solve(CHOLMOD_L, L, ZtX, &c); Free(ZtX);
    rZy = M_cholmod_solve(CHOLMOD_L, L, Zty, &c); Free(Zty);
    Memcpy(REAL(GET_SLOT(GET_SLOT(x, lme4_RZXSym), lme4_xSym)),
	   (double *) RZX->x, q * p);
    M_cholmod_free_dense(&RZX, &c);
    Memcpy(REAL(rZyP), (double *) rZy->x, q);
    M_cholmod_free_dense(&rZy, &c);
				/* signal that secondary slots are not valid */
    status[0] = 1;
}


/**
 * Downdate and factor XtX into RXX
 *
 * @param x pointer to an mer object
 *
 * @return info from the call to dpotrf
 */
static int internal_mer_Xfactor(SEXP x)
{
    int info, p = LENGTH(GET_SLOT(x, lme4_rXySym)),
	q = LENGTH(GET_SLOT(x, lme4_rZySym));
    double *RXX = REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym), lme4_xSym)),
	*RZX = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXSym), lme4_xSym)),
	one[2] = {1, 0}, m1[2] = {-1, 0};

    Memcpy(RXX, REAL(GET_SLOT(GET_SLOT(x, lme4_XtXSym), lme4_xSym)), p * p);
    F77_CALL(dsyrk)("U", "T", &p, &q, m1, RZX, &q, one, RXX, &p);
    F77_CALL(dpotrf)("U", &p, RXX, &p, &info);
    return info;
}

static void
safe_pd_matrix(double x[], const char uplo[], int n, double thresh)
{
    int info, lwork = 3 * n, nm1 = n - 1;
    double *work = Calloc(3 * n, double),
	*w = Calloc(n, double),
	*xcp = Memcpy(Calloc(n * n, double), x, n * n);

    F77_CALL(dsyev)("N", uplo, &n, xcp, &n, w, work, &lwork, &info);
    if (info) error(_("dsyev returned %d"), info);
    if (w[nm1] <= 0) error(_("no positive eigenvalues!"));
    if ((w[0]/w[nm1]) < thresh) {
	int i, np1 = n + 1;
	double incr = w[nm1] * thresh;
	for (i = 0; i < n; i++) x[i * np1] += incr;
    }
    Free(work); Free(w); Free(xcp);
}
/**
 * Update the relative precision matrices by sampling from a Wishart
 * distribution with scale factor determined by the current sample of
 * random effects.
 *
 * @param Omega pointer to the list of relative precision matrices
 * @param b current sample from the random effects
 * @param sigma current value of sigma
 * @param nf number of grouping factors
 * @param nc number columns per grouping factor
 * @param Gp integer vector pointing to the beginning of each outer
 * group of columns
 * @param vals vector in which to store values
 * @param trans logical value indicating if variance components are
 * stored in the transformed scale.
 */
static void
internal_Omega_update(SEXP Omega, const double b[], double sigma, int nf,
		  const int nc[], const int Gp[], double *vals, int trans)
{
    int i, j, k, info;
    double one = 1, zero = 0;

    for (i = 0; i < nf; i++) {
	SEXP omgi = VECTOR_ELT(Omega, i);
	int nci = nc[i];
	int nlev = (Gp[i + 1] - Gp[i])/nci, ncip1 = nci + 1,
	    ncisqr = nci * nci;
	double
	    *omgx = REAL(GET_SLOT(omgi, lme4_xSym)),
	    *scal = Calloc(ncisqr, double), /* factor of scale matrix */
	    *tmp = Calloc(ncisqr, double),
	    *var = Calloc(ncisqr, double), /* simulated variance-covariance */
	    *wfac = Calloc(ncisqr, double); /* factor of Wishart variate */
				/* create and factor the scale matrix */
	AZERO(scal, ncisqr);
	F77_CALL(dsyrk)("U", "N", &nci, &nlev, &one, b + Gp[i], &nci,
			&zero, scal, &nci);
	if (nci > 1) safe_pd_matrix(scal, "U", nci, 1e-7);
	F77_CALL(dpotrf)("U", &nci, scal, &nci, &info);
	if (info) error("Matrix is not pd after safe_pd_matrix!");
				/* generate random factor from std Wishart */
	AZERO(wfac, ncisqr);
	std_rWishart_factor((double) (nlev - nci + 1), nci, wfac);

	/* form the variance-covariance matrix and store elements */
	Memcpy(tmp, scal, ncisqr);
	F77_CALL(dtrsm)("L", "U", "T", "N", &nci, &nci,
			&one, wfac, &nci, tmp, &nci);
	F77_CALL(dsyrk)("U", "T", &nci, &nci, &one, tmp, &nci,
			&zero, var, &nci);
	for (j = 0; j < nci; j++) {
	    *vals++ = (trans ? log(var[j * ncip1]) : var[j * ncip1]);
	}
	for (j = 1; j < nci; j++) {
	    for (k = 0; k < j; k++) {
		*vals++ = (trans ? atanh(var[k + j * nci]/
				     sqrt(var[j * ncip1] * var[k * ncip1]))
		       : var[k + j * nci]);
	    }
	}
	/* calculate and store the relative precision matrix */
	Memcpy(tmp, wfac, ncisqr);
	F77_CALL(dtrsm)("R", "U", "T", "N", &nci, &nci,
			&sigma, scal, &nci, tmp, &nci);
	F77_CALL(dsyrk)("U", "T", &nci, &nci, &one, tmp, &nci, &zero,
			omgx, &nci);

	Free(scal); Free(tmp); Free(wfac); Free(var);
    }
}

/**
 * Update the fixed and random effects vectors by sampling from a
 * multivariate normal distribution with scale factor sigma, means
 * betahat and bhat, and variance-covariance matrix determined by L,
 * RZX and RXX
 *
 * @param p dimension of fixed effects vector
 * @param q dimension of random effects vector
 * @param sigma current value of sigma
 * @param L lower Cholesky factor of random-effects part of system matrix
 * @param RZX rectangular part of Cholesky factor of system matrix
 * @param RXX fixed-effectgs part of Cholesky factor of system matrix
 * @param betahat conditional mle of fixed effects
 * @param bhat conditional modes of random effects
 * @param betanew overwritten with proposal point
 * @param bnew overwritten with proposal point
 *
 * @return Gaussian kernel value for proposal point
 */
static double
internal_betab_update(int p, int q, double sigma, cholmod_factor *L,
		      double RZX[], double RXX[], double betahat[],
		      double bhat[], double betanew[], double bnew[])
{
    cholmod_dense *chb, *chbnew = M_numeric_as_chm_dense(bnew, q);
    int *perm = (int *)L->Perm;
    int j, ione = 1;
    double m1[] = {-1,0}, one[] = {1,0}, ans = 0;

				/* simulate scaled, independent normals */
    for (j = 0; j < p; j++) {
	double nr = norm_rand();
	ans += nr * nr;
	betanew[j] = sigma * nr;
    }
    for (j = 0; j < q; j++) {
	double nr = norm_rand();
	ans += nr * nr;
	bnew[perm[j]] = sigma * nr; /* RZX has permutation in it */
    }
				/* betanew := RXX^{-1} %*% betanew */
    F77_CALL(dtrsv)("U", "N", "N", &p, RXX, &p, betanew, &ione);
				/* bnew := bnew - RZX %*% betanew */
    F77_CALL(dgemv)("N", &q, &p, m1, RZX, &q, betanew, &ione,
		    one, bnew, &ione);
    for (j = 0; j < p; j++) betanew[j] += betahat[j];
				/* chb := L^{-T} %*% bnew */
    chb = M_cholmod_solve(CHOLMOD_Lt, L, chbnew, &c);
    for (j = 0; j < q; j++) {	/* Copy chb to bnew applying P-inverse */
	int pj = perm[j];
	bnew[j] = ((double*)(chb->x))[pj] + bhat[pj];
    }
    M_cholmod_free_dense(&chb, &c);
    Free(chbnew);
    return ans;
}

static double
internal_Gaussian_deviance(int p, int q, cholmod_factor *L,
			   double RZX[], double RXX[], double betahat[],
			   double bhat[], double beta[], double b[])
{
    int i, ione = 1;
    double *bb = Calloc(q, double), *betab = Calloc(p, double), ans = 0;
    cholmod_sparse *Lm;
    cholmod_dense *chb = M_numeric_as_chm_dense(bb, q),
	*Ltb = M_cholmod_allocate_dense(q, 1, q, CHOLMOD_REAL, &c);
    cholmod_factor *Lcp;
    double one[] = {1,0}, zero[] = {0,0};

    for (i = 0; i < p; i++) betab[i] = beta[i] - betahat[i];
    for (i = 0; i < q; i++) bb[i] = b[i] - bhat[i];
    Lcp = M_cholmod_copy_factor(L, &c); /* next call changes Lcp */
    Lm = M_cholmod_factor_to_sparse(Lcp, &c); M_cholmod_free_factor(&Lcp, &c);
    if (!M_cholmod_sdmult(Lm, 1 /* transpose */, one, zero, chb, Ltb, &c))
	error(_("Error return from cholmod_sdmult"));
    Memcpy(bb, (double *)(Ltb->x), q);
    M_cholmod_free_sparse(&Lm, &c); M_cholmod_free_dense(&Ltb, &c);
    F77_CALL(dgemv)("N", &q, &p, one, RZX, &q, betab, &ione, one, bb, &ione);
    for (i = 0; i < q; i++) ans += bb[i] * bb[i];

    F77_CALL(dtrmv)("U", "N", "N", &p, RXX, &p, betab, &ione);
    for (i = 0; i < p; i++) ans += betab[i] * betab[i];
    Free(chb); Free(bb); Free(betab);
    return ans;
}


/**
 * Evaluate new weights and working residuals.
 *
 * @param GS a GlmerStruct object
 */
static void
internal_glmer_reweight(GlmerStruct GS) {
    SEXP dmu_deta, var;
    int i;
    double *w = REAL(GET_SLOT(GS->mer, lme4_wtsSym)),
	*y = REAL(GET_SLOT(GS->mer, lme4_ySym)),
	*z = REAL(GET_SLOT(GS->mer, lme4_wrkresSym));

				/* reweight mer */
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    dmu_deta = PROTECT(eval_check(GS->mu_eta, GS->rho,
				  REALSXP, GS->n));
    var = PROTECT(eval_check(GS->var, GS->rho,
			     REALSXP, GS->n));
    /* FIXME: speedup doing REAL(.) outside of loop :*/
    for (i = 0; i < GS->n; i++) {
	w[i] = GS->wts[i] * REAL(dmu_deta)[i]/sqrt(REAL(var)[i]);
	z[i] = REAL(GS->eta)[i] - GS->offset[i] +
	    (y[i] - REAL(GS->mu)[i])/REAL(dmu_deta)[i];
    }
    UNPROTECT(2);
    mer_update_ZXy(GS->mer);
}

/**
 * Update eta, evaluate the convergence criterion, then copy eta to
 * etaold
 *
 * @param GS a GlmerStruct object
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
 * Update the ranef slot assuming that the fixef, rZy, RZX and L slots
 * are up to date.
 *
 * @param x Pointer to an mer object
 *
 * @return
 */
static double*
internal_mer_ranef(SEXP x)
{
    SEXP ranef = GET_SLOT(x, lme4_ranefSym);
    int *status = INTEGER(GET_SLOT(x, lme4_statusSym));
    if (!status[0]) {
	error("Applying internal_mer_ranef without factoring");
	return (double*)NULL;	/* -Wall */
    }
    if (status[0] < 2) {
	SEXP fixef = GET_SLOT(x, lme4_fixefSym),
	    ranef = GET_SLOT(x, lme4_ranefSym);
	int ione = 1, p = LENGTH(fixef), q = LENGTH(ranef);
	cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
	cholmod_dense *td1, *td2,
	    *chranef = M_as_cholmod_dense(ranef);
	double *RZX = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXSym), lme4_xSym)),
	    m1[] = {-1,0}, one[] = {1,0};

	Memcpy(REAL(ranef), REAL(GET_SLOT(x, lme4_rZySym)), q);
	F77_CALL(dgemv)("N", &q, &p, m1, RZX, &q, REAL(fixef), &ione,
			one, REAL(ranef), &ione);
	td1 = M_cholmod_solve(CHOLMOD_Lt, L, chranef, &c);
	td2 = M_cholmod_solve(CHOLMOD_Pt, L, td1, &c);
	Free(chranef); M_cholmod_free_dense(&td1, &c);
	Memcpy(REAL(ranef), (double *)(td2->x), q);
	M_cholmod_free_dense(&td2, &c);
	status[0] = 2;
	Free(L);
    }
    return REAL(ranef);
}

/**
 * Update the fixef slot on a factored mer object.
 *
 * @param x Pointer to an mer object
 *
 * @return fixed effects vector
 */
static double*
internal_mer_fixef(SEXP x)
{
    SEXP fixef = GET_SLOT(x, lme4_fixefSym);
    int *status = INTEGER(GET_SLOT(x, lme4_statusSym));
    if (!status[0]) {
	error("Applying internal_mer_fixef without factoring");
	return (double*)NULL;	/* -Wall */
    }
    if (status[0] < 2) {
	int ione = 1, p = LENGTH(fixef);
	Memcpy(REAL(fixef), REAL(GET_SLOT(x, lme4_rXySym)), p);
	F77_CALL(dtrsv)("U", "N", "N", &p,
			REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym),
				      lme4_xSym)),
			&p, REAL(fixef), &ione);
    }
    return REAL(fixef);
}

/**
 * Iterate to determine the conditional modes of the random effects.
 *
 * @param GS a GlmerStruct object
 * @param fixed vector of fixed effects
 * @param varc vector of parameters for the variance-covariance
 *
 * @return An indicator of whether the iterations converged
 */
static int
internal_bhat(GlmerStruct GS, const double fixed[], const double varc[])
{
    SEXP fixef = GET_SLOT(GS->mer, lme4_fixefSym);
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(GS->mer, lme4_LSym));
    int i;
    double crit = GS->tol + 1, *etap = REAL(GS->eta);

    if (varc)	  /* skip this step if varc == (double*) NULL */
	internal_mer_coefGets(GS->mer, varc, 2);
    Memcpy(REAL(fixef), fixed, LENGTH(fixef));
    internal_mer_Zfactor(GS->mer, L);
    internal_mer_ranef(GS->mer);
    internal_mer_fitted(GS->mer, GS->offset, etap);
    Memcpy(GS->etaold, etap, GS->n);

    for (i = 0; i < GS->maxiter && crit > GS->tol; i++) {
	internal_glmer_reweight(GS);
	internal_mer_Zfactor(GS->mer, L);
	internal_mer_ranef(GS->mer);
	internal_mer_fitted(GS->mer, GS->offset, etap);
	crit = conv_crit(GS->etaold, etap, GS->n);
    }
    internal_mer_Xfactor(GS->mer);
    Free(L);
    return (crit > GS->tol) ? 0 : i;
}

/**
 * Print the verbose output in the ECME iterations
 *
 * @param x pointer to an ssclme object
 * @param iter iteration number
 * @param REML non-zero for REML, zero for ML
 */
static void
EMsteps_verbose_print(SEXP x, int iter, int REML)
{
    SEXP Omega = GET_SLOT(x, lme4_OmegaSym),
	gradComp = GET_SLOT(x, lme4_gradCompSym);
    int *nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	i, ifour = 4, ii, ione = 1, jj,
	nf = LENGTH(Omega);
    double
	cc[] = {-1, 1, 1, REML ? 1 : 0},
	*dev = REAL(GET_SLOT(x, lme4_devianceSym)),
	one = 1., zero = 0.;

    if (iter == 0) Rprintf("  EM iterations\n");
    Rprintf("%3d %.3f", iter, dev[REML ? 1 : 0]);
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1, ncisqr = nci * nci;
	double
	    *Omgi = REAL(GET_SLOT(VECTOR_ELT(Omega, i), lme4_xSym)),
	    *Grad = Calloc(ncisqr, double);

				/* diagonals */
	Rprintf(" (%#8g", Omgi[0]);
	for (jj = 1; jj < nci; jj++) {
	    Rprintf(" %#8g", Omgi[jj * ncip1]);
	}
	for (jj = 1; jj < nci; jj++) /* offdiagonals */
	    for (ii = 0; ii < jj; ii++)
		Rprintf(" %#8g", Omgi[ii + jj * nci]);
				/* Evaluate and print the gradient */
	F77_CALL(dgemv)("N", &ncisqr, &ifour, &one,
			REAL(VECTOR_ELT(gradComp, i)), &ncisqr,
			cc, &ione, &zero, Grad, &ione);
	Rprintf(":%#8.3g", Grad[0]);

	for (jj = 1; jj < nci; jj++) { /* diagonals */
	    Rprintf(" %#8.3g", Grad[jj * ncip1]);
	}
	for (jj = 1; jj < nci; jj++) /* offdiagonals */
	    for (ii = 0; ii < jj; ii++)
		Rprintf(" %#8.3g", Grad[ii + jj * nci]);
	Rprintf(")");
	Free(Grad);
    }
    Rprintf("\n");
}

/**
 * Perform a number of ECME steps
 *
 * @param x pointer to an mer object
 * @param nEM number of iteration to perform
 * @param verb indicator of verbose output
 */
static void
internal_ECMEsteps(SEXP x, int nEM, int verb)
{
    SEXP Omega = GET_SLOT(x, lme4_OmegaSym),
	gradComp = GET_SLOT(x, lme4_gradCompSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	REML = INTEGER(GET_SLOT(x, lme4_statusSym))[1],
	i, ifour = 4, info, ione = 1, iter,
	nf = LENGTH(Omega);
    double
	cc[] = {0, 1, 1, REML ? 1 : 0},
	zero = 0.0;

    mer_gradComp(x);
    if (verb)
	EMsteps_verbose_print(x, 0, REML);
    for (iter = 0; iter < nEM; iter++) {
	for (i = 0; i < nf; i++) {
	    SEXP Omgi = VECTOR_ELT(Omega, i);
	    int nci = nc[i], ncisqr = nci * nci,
		nlev = (Gp[i + 1] - Gp[i])/nci;
	    double *Omgx = REAL(GET_SLOT(Omgi, lme4_xSym)),
		mult = 1./((double) nlev);

	    F77_CALL(dgemm)("N", "N", &ncisqr, &ione, &ifour, &mult,
			    REAL(VECTOR_ELT(gradComp, i)), &ncisqr,
			    cc, &ifour, &zero, Omgx, &ncisqr);
	    F77_CALL(dpotrf)("U", &nci, Omgx, &nci, &info);
	    if (info)
		error(_("DPOTRF in ECME update gave code %d"),
		      info);
	    F77_CALL(dpotri)("U", &nci, Omgx, &nci, &info);
	    if (info)
		error(_("Matrix inverse in ECME update gave code %d"), info);
	    SET_SLOT(Omgi, lme4_factorSym, allocVector(VECSXP, 0));
	}
	flag_not_factored(x);
	mer_gradComp(x);
	if (verb)
	    EMsteps_verbose_print(x, iter + 1, REML);
    }
}

/**
 * Evaluate the conditional deviance for the stored random effects.
 *
 * @param GS Pointer to a GlmerStruct
 *
 * @return conditional deviance
 */
static double
random_effects_deviance(GlmerStruct GS)
{
    int i;
    double ans, *devs;

    internal_mer_fitted(GS->mer, GS->offset, REAL(GS->eta));
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    devs = REAL(PROTECT(eval_check(GS->dev_resids, GS->rho, REALSXP, GS->n)));
    for (i = 0, ans = 0; i < GS->n; i++) ans += devs[i];
    UNPROTECT(1);
    return ans;
}

/**
 * Calculate the deviance for a generalizedlinear mixed model at
 * arbitrary parameter values.  This version restores the original
 * values of the fixef and ranef slots after evaluating at arbitrary
 * beta and b.
 *
 * @param GS a generalized mixed-effects model pointer
 * @param beta fixed-effects parameter vector
 * @param b random-effects vector
 *
 * @return deviance
 */
static
double glmm_deviance(GlmerStruct GS, const double beta[], const double b[])
{
    SEXP x = GS->mer;
    SEXP fixefp = GET_SLOT(x, lme4_fixefSym),
	ranefp = GET_SLOT(x, lme4_ranefSym);
    int p = LENGTH(fixefp), q = LENGTH(ranefp);
    double *fixcp = Memcpy(Calloc(p, double), REAL(fixefp), p),
	*rancp = Memcpy(Calloc(q, double), REAL(ranefp), q), ans;

    mer_factor(x);
    Memcpy(REAL(fixefp), beta, p);
    Memcpy(REAL(ranefp), b, q);
    ans = random_effects_deviance(GS) +
	b_quadratic(b, GET_SLOT(x, lme4_OmegaSym),
		    INTEGER(GET_SLOT(x, lme4_GpSym)),
		    INTEGER(GET_SLOT(x, lme4_ncSym)));
    Memcpy(REAL(fixefp), fixcp, p);
    Memcpy(REAL(ranefp), rancp, q);
    Free(fixcp); Free(rancp);
    return ans;
}


/* Externally accessible functions */

/**
 * Simulate a sample of random matrices from a Wishart distribution
 *
 * @param ns Number of samples to generate
 * @param dfp Degrees of freedom
 * @param scal Positive-definite scale matrix
 *
 * @return
 */
SEXP
lme4_rWishart(SEXP ns, SEXP dfp, SEXP scal)
{
    SEXP ans;
    int *dims = INTEGER(getAttrib(scal, R_DimSymbol)), j,
	n = asInteger(ns), psqr;
    double *scCp, *ansp, *tmp, df = asReal(dfp), one = 1, zero = 0;

    if (!isMatrix(scal) || !isReal(scal) || dims[0] != dims[1])
	error("scal must be a square, real matrix");
    if (n <= 0) n = 1;
    psqr = dims[0] * dims[0];
    tmp = Calloc(psqr, double);
    AZERO(tmp, psqr);
    scCp = Memcpy(Calloc(psqr, double), REAL(scal), psqr);
    F77_CALL(dpotrf)("U", &(dims[0]), scCp, &(dims[0]), &j);
    if (j)
	error("scal matrix is not positive-definite");
    PROTECT(ans = alloc3Darray(REALSXP, dims[0], dims[0], n));
    ansp = REAL(ans);
    GetRNGstate();
    for (j = 0; j < n; j++) {
	double *ansj = ansp + j * psqr;
	std_rWishart_factor(df, dims[0], tmp);
	F77_CALL(dtrmm)("R", "U", "N", "N", dims, dims,
			&one, scCp, dims, tmp, dims);
	F77_CALL(dsyrk)("U", "T", &(dims[1]), &(dims[1]),
			&one, tmp, &(dims[1]),
			&zero, ansj, &(dims[1]));
	internal_symmetrize(ansj, dims[0]);
    }

    PutRNGstate();
    Free(scCp); Free(tmp);
    UNPROTECT(1);
    return ans;
}

/**
 * Perform the PQL optimization
 *
 * @param GSp pointer to a GlmerStruct object
 *
 * @return R_NilValue
 */
SEXP glmer_PQL(SEXP GSp)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    int i; double crit, *etap = REAL(GS->eta);

    Memcpy(GS->etaold, etap, GS->n);
    for (i = 0, crit = GS->tol + 1;
	 i < GS->maxiter && crit > GS->tol; i++) {
	internal_glmer_reweight(GS);
	if (!i) mer_initial(GS->mer); /* initialize first fit */
	internal_ECMEsteps(GS->mer, i ? 2 : GS->niterEM,
			   GS->EMverbose);
	eval(GS->LMEopt, GS->rho);
	internal_mer_fitted(GS->mer, GS->offset, etap);
	crit = conv_crit(GS->etaold, etap, GS->n);
    }
    if (crit > GS->tol)
	warning(_("IRLS iterations for PQL did not converge"));

    return R_NilValue;
}

/**
 * Compute the Laplace approximation to the deviance.
 *
 * @param pars pointer to a numeric vector of parameters
 * @param GSp pointer to a GlmerStruct object
 *
 * @return the Laplace approximation to the deviance
 */
SEXP glmer_devLaplace(SEXP pars, SEXP GSp)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    SEXP Omega = GET_SLOT(GS->mer, lme4_OmegaSym);
    int *Gp = INTEGER(GET_SLOT(GS->mer, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(GS->mer, lme4_ncSym));
    double *bhat = REAL(GET_SLOT(GS->mer, lme4_ranefSym)),
	*dcmp = REAL(GET_SLOT(GS->mer, lme4_devCompSym)),
	*dev = REAL(GET_SLOT(GS->mer, lme4_devianceSym));

    if (!isReal(pars) || LENGTH(pars) != GS->npar)
	error(_("`%s' must be a numeric vector of length %d"),
	      "pars", GS->npar);
    if (!internal_bhat(GS, REAL(pars), REAL(pars) + (GS->p)))
	return ScalarReal(DBL_MAX);
    dev[0] = dcmp[4] - dcmp[5] + random_effects_deviance(GS) +
	b_quadratic(bhat, Omega, Gp, nc);
    dev[1] = NA_REAL;
    return ScalarReal(dev[0]);
}

/**
 * Release the storage for a GlmerStruct
 *
 * @param GSp External pointer to a  GlmerStruct
 *
 * @return R_NilValue
 */
SEXP glmer_finalize(SEXP GSp) {
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);

    Free(GS->offset); Free(GS->wts); Free(GS->etaold);
    Free(GS);
    return R_NilValue;
}

/**
 * Return an external pointer object to a GlmerStruct created in
 * environment rho
 *
 * @param rho An environment
 *
 * @return An external pointer to a GlmerStruct
 */
SEXP glmer_init(SEXP rho) {
    GlmerStruct GS;
    SEXP tmp, y, Ztx;


    GS = (GlmerStruct) Calloc(1, glmer_struct);
    if (!isEnvironment(rho))
	error(_("`rho' must be an environment"));
    GS->rho = rho;
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 4, 0)
    GS->mer = find_and_check(rho, install("mer"), S4SXP, 0);
#else
    GS->mer = find_and_check(rho, install("mer"), VECSXP, 0);
#endif /* S4SXP */
    y = GET_SLOT(GS->mer, lme4_ySym);
    GS->n = LENGTH(y);
    GS->p = LENGTH(GET_SLOT(GS->mer, lme4_rXySym));
    GS->y = Memcpy(Calloc(GS->n, double), REAL(y), GS->n);
    Ztx = GET_SLOT(GET_SLOT(GS->mer, lme4_ZtSym), lme4_xSym);
    GS->eta = find_and_check(rho, install("eta"), REALSXP, GS->n);
    GS->mu = find_and_check(rho, install("mu"), REALSXP, GS->n);
    tmp = find_and_check(rho, install("offset"), REALSXP, GS->n);
    GS->offset = Memcpy(Calloc(GS->n, double), REAL(tmp), GS->n);
    tmp = find_and_check(rho, install("weights"), REALSXP, GS->n);
    GS->wts = Memcpy(Calloc(GS->n, double), REAL(tmp), GS->n);
    GS->etaold = Calloc(GS->n, double);
    GS->cv = find_and_check(rho, install("cv"), VECSXP, 0);
    GS->niterEM = asInteger(internal_getElement(GS->cv, "niterEM"));
    GS->EMverbose = asLogical(internal_getElement(GS->cv, "EMverbose"));
    GS->tol = asReal(internal_getElement(GS->cv, "tolerance"));
    GS->maxiter = asInteger(internal_getElement(GS->cv, "maxIter"));
    GS->nf = LENGTH(GET_SLOT(GS->mer, lme4_flistSym));
    GS->npar = GS->p +
	coef_length(GS->nf, INTEGER(GET_SLOT(GS->mer, lme4_ncSym)));

    GS->linkinv = find_and_check(rho, install("linkinv"),
				 LANGSXP, 0);
    GS->mu_eta = find_and_check(rho, install("mu.eta"),
				LANGSXP, 0);
    GS->var = find_and_check(rho, install("variance"),
			     LANGSXP, 0);
    GS->LMEopt = find_and_check(rho, install("doLMEopt"),
				LANGSXP, 0);
    GS->dev_resids = find_and_check(rho, install("dev.resids"),
				    LANGSXP, 0);
    return R_MakeExternalPtr(GS, R_NilValue, GS->mer);
}


/**
 * Perform ECME steps for the REML or ML criterion.
 *
 * @param x pointer to an mer object
 * @param nsteps pointer to an integer scalar - the number of ECME
 * steps to perform
 * @param Verbp pointer to a logical scalar indicating verbose output
 *
 * @return R_NilValue
 */
SEXP mer_ECMEsteps(SEXP x, SEXP nsteps, SEXP Verbp)
{
    int nstp = asInteger(nsteps);
    if (nstp > 0) internal_ECMEsteps(x, nstp, asLogical(Verbp));
    return R_NilValue;
}

/**
 * Fill in five symmetric matrices, providing the information to
 * generate the Hessian.
 *
 * @param x pointer to an mer object
 *
 * @return an array consisting of five symmetric faces
 */
SEXP mer_Hessian(SEXP x)
{
    SEXP
	D = GET_SLOT(x, lme4_DSym),
	Omega = GET_SLOT(x, lme4_OmegaSym),
	RZXP = GET_SLOT(x, lme4_RZXSym),
	gradComp = GET_SLOT(x, lme4_gradCompSym),
	val;
    int *dRZX = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	*Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	Q, Qsqr, RZXpos, facepos,
	i, ione = 1, j, nf = length(Omega), p = dRZX[1] - 1, pos;
    double
	*RZX = REAL(RZXP),
	*b = REAL(RZXP) + dRZX[0] * p,
	*valp, *bbface,		/* vec of second faces of gradComp elts */
	one = 1.,
	zero = 0.;

    mer_gradComp(x);
    Q = 0;			/* number of rows and columns in the result */
    for (i = 0; i < nf; i++) Q += nc[i] * nc[i];
    Qsqr = Q * Q;
    bbface = Calloc(Q, double);
    val = PROTECT(alloc3Darray(REALSXP, Q, Q, 5));
    valp= REAL(val);
    AZERO(valp, Qsqr * 5);

    pos = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	int ncisqr = nci * nci;
	double *fDi = REAL(VECTOR_ELT(gradComp, i)),
	    mult = 1./((double)(Gp[i + 1] - Gp[i])/nci);

	Memcpy(bbface + pos, fDi + ncisqr, ncisqr);
	/* outer product of the third face of gradComp on the diagonal
	 * of the third face of val */
	F77_CALL(dsyr)("U", &ncisqr, &mult, fDi + 2 * ncisqr, &ione,
		       valp + 2 * Qsqr + pos * Q, &Q);
	pos += ncisqr;
    }
				/* fifth face is outer product of bbface */
    F77_CALL(dsyr)("U", &Q, &one, bbface, &ione, valp + 4 * Qsqr, &Q);
				/* fourth face from \bb\trans\der\vb\der\bb */
    AZERO(valp + 3 * Qsqr, Qsqr); /* zero accumulator */
    RZXpos = 0;
    facepos = 0;
    for (i = 0; i < nf; i++) {
	int ii, jj, nci = nc[i];
	int ncisqr = nci * nci, nctp = nci * p,
	    nlev = (Gp[i + 1] - Gp[i])/nci;
	int maxpq = (p > nci) ? p : nci;
	double
	    *Di = REAL(VECTOR_ELT(D, i)),
	    *arr = Calloc(ncisqr * maxpq, double), /* tmp 3Darray */
	    *face = valp + 3 * Qsqr,
	    *mat = Calloc(nci * maxpq, double); /* tmp matrix */

	for (j = 0; j < nlev; j++) {
	    F77_CALL(dgemm)("T", "T", &p, &nci, &nci,
			    &one, RZX + j * nci, dRZX, Di + j * ncisqr, &nci,
			    &zero, mat, &p);
	    F77_CALL(dgemm)("N", "N", &nctp, &nci, &ione,
			    &one, mat, &nctp, b + j * nci, &ione,
			    &zero, arr, &nctp);
	    F77_CALL(dsyrk)("U", "T", &ncisqr, &p, &one, arr, &p,
			    &one, face + facepos, &Q);
				/* Add the D_{i,j}^{-T/2} term */
	    Memcpy(mat, Di + j * ncisqr, ncisqr);
	    for (jj = 1; jj < nci; jj++) { /* transpose mat */
		for (ii = 0; ii < jj; ii++) {
		    mat[jj + ii * nci] = mat[ii + jj * nci];
		    mat[ii + jj * nci] = 0.;
		}
	    }
	    F77_CALL(dgemm)("N", "N", &ncisqr, &nci, &ione,
			    &one, mat, &ncisqr, b + j * nci, &ione,
			    &zero, arr, &ncisqr);
	    /* FIXME: Next call could be dsyr (it's rank one). */
	    F77_CALL(dsyrk)("U", "T", &ncisqr, &nci, &one, arr, &nci,
			    &one, face + facepos, &Q);

	}
	RZXpos += nci * nlev;
	facepos += ncisqr;
	Free(arr); Free(mat);
    }
    UNPROTECT(2);
    Free(bbface);
    return val;
}

/**
 * Generate a Markov-Chain Monte Carlo sample from a fitted
 * linear mixed model.
 *
 * @param x pointer to an mer object
 * @param savebp pointer to a logical scalar indicating if the
 * random-effects should be saved
 * @param nsampp pointer to an integer scalar of the number of samples
 * to generate
 * @param transp pointer to an logical scalar indicating if the
 * variance components should be transformed.
 *
 * @return a matrix
 */
SEXP mer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp, SEXP verbosep)
{
    SEXP ans, Omega = GET_SLOT(x, lme4_OmegaSym),
	Omegacp = PROTECT(duplicate(Omega));
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	REML = INTEGER(GET_SLOT(x, lme4_statusSym))[1],
	i, j, n = LENGTH(GET_SLOT(x, lme4_ySym)),
	nf = LENGTH(Omega), nsamp = asInteger(nsampp),
	p = LENGTH(GET_SLOT(x, lme4_rXySym)),
	q = LENGTH(GET_SLOT(x, lme4_rZySym)),
	saveb = asLogical(savebp),
	trans = asLogical(transp),
	verbose = asLogical(verbosep);/* currently unused */
    double
	*RXX = REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym), lme4_xSym)),
	*RZX = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXSym), lme4_xSym)),
	*bhat = REAL(GET_SLOT(x, lme4_ranefSym)),
	*betahat = REAL(GET_SLOT(x, lme4_fixefSym)),
	*bnew = Calloc(q, double), *betanew = Calloc(p, double),
	*dcmp = REAL(GET_SLOT(x, lme4_devCompSym)),
	*ansp, df = n - (REML ? p : 0);
    int nrbase = p + 2 + coef_length(nf, nc); /* rows always included */
    int nrtot = nrbase + (saveb ? q : 0);
    cholmod_dense *chbnew = M_numeric_as_chm_dense(bnew, q);

    if (nsamp <= 0) nsamp = 1;
    ans = PROTECT(allocMatrix(REALSXP, nrtot, nsamp));
    ansp = REAL(ans);
    for (i = 0; i < nrtot * nsamp; i++) ansp[i] = NA_REAL;
    GetRNGstate();
    if(verbose) Rprintf("%12s %14s\n", "sigma", "deviance");

    for (i = 0; i < nsamp; i++) {
	double *col = ansp + i * nrtot, sigma, dev;
				/* simulate and store new value of sigma */
	sigma = exp(dcmp[3]/2)/sqrt(rchisq(df));
	col[p] = (trans ? 2 * log(sigma) : sigma * sigma);
				/* simulate new fixed and random effects */
	internal_betab_update(p, q, sigma, L, RZX, RXX, betahat, bhat,
			      betanew, bnew);
				/* Store beta */
	for (j = 0; j < p; j++) col[j] = betanew[j];
				/* Optionally store b */
	if (saveb) for (j = 0; j < q; j++) col[nrbase + j] = bnew[j];
				/* Update and store Omega */
	internal_Omega_update(Omega, bnew, sigma, nf, nc, Gp,
				     col + p + 1, trans);
	internal_mer_refactor(x);
	mer_secondary(x);

	dev = lmm_deviance(x, sigma, betanew);
	col[nrbase - 1] = dev; /* store deviance */
	if(verbose) Rprintf("%12.6g %14.8g\n", sigma, dev);
    }
    PutRNGstate();
    Free(betanew); Free(bnew); Free(chbnew);
				/* Restore original Omega */
    SET_SLOT(x, lme4_OmegaSym, Omegacp);
    internal_mer_refactor(x);

    Free(L);
    UNPROTECT(2);
    return ans;
}

/**
 * Create an mer object from a list of grouping factors and a list of model
 * matrices.
 *
 * @param fl named list of grouping factors
 * @param ZZt transpose of Z as a sparse matrix
 * @param Xp model matrix for the fixed effects
 * @param yp response vector
 * @param REMLp logical scalar indicating if REML is to be used
 * @param ncp integer vector of the number of random effects per level
 *        of each grouping factors
 *
 * @return pointer to an mer object
 */
SEXP mer_create(SEXP fl, SEXP ZZt, SEXP Xp, SEXP yp, SEXP REMLp,
		SEXP ncp, SEXP cnames)
{
    SEXP Omega, bVar, gradComp, fnms = getAttrib(fl, R_NamesSymbol),
	stat, val = PROTECT(NEW_OBJECT(MAKE_CLASS("mer"))), xnms;
    cholmod_sparse *ts1, *ts2, *Zt;
    cholmod_factor *F;
    int *nc = INTEGER(ncp), *Gp, *xdims, REML = asInteger(REMLp),
	i, nested, nf = LENGTH(fl), nobs = LENGTH(yp), p, q;
    char *devCmpnms[] = {"n", "p", "yty", "logryy2", "logDetL2",
			 "logDetOmega", "logDetRXX", "scale", ""};
    char *devnms[] = {"ML", "REML", ""};
    char *statusnms[] =
	{"stage", "REML", "glmm", ""};
    double *dcmp, *wts, *wrkres, *y;
				/* Check arguments to be duplicated */
    if (!isReal(yp)) error(_("yp must be a real vector"));
    SET_SLOT(val, lme4_ySym, duplicate(yp));
    if (!isMatrix(Xp) || !isReal(Xp))
	error(_("Xp must be a real matrix"));
    xdims = INTEGER(getAttrib(Xp, R_DimSymbol));
    if (xdims[0] != nobs) error(_("Xp must have %d rows"), nobs);
    p = xdims[1];
    xnms = VECTOR_ELT(getAttrib(Xp, R_DimNamesSymbol), 1);
    SET_SLOT(val, lme4_XSym, duplicate(Xp));
    if (!isNewList(fl) || nf < 1) error(_("fl must be a nonempty list"));
    for (i = 0; i < nf; i++) {
	SEXP fli = VECTOR_ELT(fl, i);
	if (!isFactor(fli) || LENGTH(fli) != nobs)
	    error(_("fl[[%d] must be a factor of length %d"), i+1, nobs);
    }
    SET_SLOT(val, lme4_flistSym, duplicate(fl));
    if (!isNewList(cnames) || LENGTH(cnames) != nf + 1)
	error(_("cnames must be a list of length %d"), nf + 1);
    SET_SLOT(val, lme4_cnamesSym, duplicate(cnames));
    if (!isInteger(ncp) || LENGTH(ncp) != nf)
	error(_("ncp must be an integer vector of length %d"), nf);
    SET_SLOT(val, lme4_ncSym, duplicate(ncp));
    Gp = INTEGER(ALLOC_SLOT(val, lme4_GpSym, INTSXP, nf + 1));
    Gp[0] = 0;
    if (!isNewList(fl) || nf < 1) error(_("fl must be a nonempty list"));
    for (i = 0; i < nf; i++) {
	SEXP fli = VECTOR_ELT(fl, i);
	if (!isFactor(fli) || LENGTH(fli) != nobs)
	    error(_("fl[[%d] must be a factor of length %d"), i+1, nobs);

    }
    SET_SLOT(val, lme4_ZtSym, duplicate(ZZt));
    Zt = M_as_cholmod_sparse(GET_SLOT(val, lme4_ZtSym));
    q = Zt->nrow;
				/* allocate other slots */
    SET_SLOT(val, lme4_devianceSym, internal_make_named(REALSXP, devnms));
    SET_SLOT(val, lme4_devCompSym, internal_make_named(REALSXP, devCmpnms));
    SET_SLOT(val, lme4_statusSym, internal_make_named(INTSXP, statusnms));
    stat = GET_SLOT(val, lme4_statusSym);
    AZERO(INTEGER(stat), LENGTH(stat));
    INTEGER(stat)[1] = REML;
    dcmp = REAL(GET_SLOT(val, lme4_devCompSym));
    AZERO(dcmp, 8);		/* cosmetic */
    dcmp[0] = (double) nobs;
    dcmp[1] = (double) p;
    dcmp[7] = 1.;		/* initialize to a positive value */
				/* allocate and populate list slots */
    Omega = ALLOC_SLOT(val, lme4_OmegaSym, VECSXP, nf);
    bVar = ALLOC_SLOT(val, lme4_bVarSym, VECSXP, nf);
    gradComp = ALLOC_SLOT(val, lme4_gradCompSym, VECSXP, nf);
    setAttrib(Omega, R_NamesSymbol, duplicate(fnms));
    setAttrib(bVar, R_NamesSymbol, duplicate(fnms));
    setAttrib(gradComp, R_NamesSymbol, duplicate(fnms));
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	int nlev = LENGTH(getAttrib(VECTOR_ELT(fl, i), R_LevelsSymbol));
	SET_VECTOR_ELT(Omega, i,
		       alloc_dpoMatrix(nci, "U",
					 VECTOR_ELT(cnames, i),
					 VECTOR_ELT(cnames, i)));
	SET_VECTOR_ELT(bVar, i, alloc3Darray(REALSXP, nci, nci, nlev));
	SET_VECTOR_ELT(gradComp, i, alloc3Darray(REALSXP, nci, nci, 4));
	Gp[i + 1] = Gp[i] + nlev * nci;
    }
				/* analyze Zt and ZtZ */
    ts1 = M_cholmod_aat(Zt, (int*)NULL/* fset */,(size_t)0,
		      CHOLMOD_PATTERN, &c);
    ts2 = M_cholmod_copy(ts1, -1/*lower triangle*/, CHOLMOD_PATTERN, &c);
    SET_SLOT(val, lme4_ZtZSym,
	     alloc_dsCMatrix(q, M_cholmod_nnz(ts2, &c), "U", R_NilValue,
			       R_NilValue));
    i = c.supernodal;
    c.supernodal = CHOLMOD_SUPERNODAL; /* force a supernodal decomposition */
    nested = internal_mer_isNested(nf, nc, Gp, (int*)(ts2->p));
    if (nested) {		/* require identity permutation */
	int nmethods = c.nmethods, ord0 = c.method[0].ordering;
	c.nmethods = 1;
	c.method[0].ordering = CHOLMOD_NATURAL;
	c.postorder = FALSE;
	F = M_cholmod_analyze(Zt, &c);
	c.nmethods = nmethods; c.method[0].ordering = ord0;
	c.postorder = TRUE;
    } else F = M_cholmod_analyze(Zt, &c);
    c.supernodal = i;		/* restore previous setting */
    M_cholmod_free_sparse(&ts1, &c); M_cholmod_free_sparse(&ts2, &c);
				/* create ZtX, RZX, XtX, RXX */
    SET_SLOT(val, lme4_ZtXSym, alloc_dgeMatrix(q, p, R_NilValue, xnms));
    SET_SLOT(val, lme4_RZXSym, alloc_dgeMatrix(q, p, R_NilValue, xnms));
    SET_SLOT(val, lme4_XtXSym, alloc_dpoMatrix(p, "U", xnms, xnms));
    SET_SLOT(val, lme4_RXXSym, alloc_dtrMatrix(p, "U", "N", xnms, xnms));
    SET_SLOT(val, lme4_ZtySym, allocVector(REALSXP, q));
    SET_SLOT(val, lme4_rZySym, allocVector(REALSXP, q));
    SET_SLOT(val, lme4_XtySym, allocVector(REALSXP, p));
    SET_SLOT(val, lme4_rXySym, allocVector(REALSXP, p));
				/* create weights and working residuals */
    wts = REAL(ALLOC_SLOT(val, lme4_wtsSym, REALSXP, nobs));
    wrkres = REAL(ALLOC_SLOT(val, lme4_wrkresSym, REALSXP, nobs));
    y = REAL(GET_SLOT(val, lme4_ySym));
    for (i = 0; i < nobs; i++) {wts[i] = 1.; wrkres[i] = y[i];}
    internal_mer_update_ZXy(val, (int*)(F->Perm));
    Free(Zt);
				/* secondary slots */
    SET_SLOT(val, lme4_ranefSym, allocVector(REALSXP, q));
    SET_SLOT(val, lme4_fixefSym, allocVector(REALSXP, p));
    SET_SLOT(val, lme4_RZXinvSym, alloc_dgeMatrix(q, p, R_NilValue, xnms));
				/* initialize */
    mer_initial(val);
    /* The next calls are simply to set up the L slot.  At present the
     * factor F is a symbolic factor.  We need to force it to become
     * numeric before allocating the L slot in the object. */
    internal_mer_Zfactor(val, F);
    /* One side-effect of this call is to set the status as
     * factored.  We need to reset it */
    INTEGER(GET_SLOT(val, lme4_statusSym))[0] = 0;
    /* Create the dCHMfactor object and store it in the L slot.  This
     * also frees the storage. */
    SET_SLOT(val, lme4_LSym, M_chm_factor_to_SEXP(F, 1));
    /* OK, done now. */
    UNPROTECT(1);
    return val;
}

/**
 * Extract parameters from the Omega matrices.  These aren't
 * "coefficients" but the extractor is called coef for historical
 * reasons.  Within each group these values are in the order of the
 * diagonal entries first then the strict upper triangle in row
 * order.
 *
 * The parameters can be returned in three forms:
 *   0 - nonlinearly constrained - elements of the relative precision matrix
 *   1 - unconstrained - from the LDL' decomposition - logarithms of
 *       the diagonal elements of D
 *   2 - box constrained - also from the LDL' decomposition - inverses
 *       of the diagonal elements of D
 *
 * @param x pointer to an mer object
 * @param pType pointer to an integer scalar indicating the form of the
 *        parameters to be returned.
 *
 * @return numeric vector of the values in the upper triangles of the
 * Omega matrices
 */
SEXP mer_coef(SEXP x, SEXP pType)
{
    int	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	nf = LENGTH(GET_SLOT(x, lme4_OmegaSym));
    SEXP val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));

    internal_mer_coef(x, asInteger(pType), REAL(val));
    UNPROTECT(1);
    return val;
}

/**
 * Assign the upper triangles of the Omega matrices according to a
 * vector of parameters.
 *
 * @param x pointer to an lme object
 * @param coef pointer to an numeric vector of appropriate length
 * @param pType pointer to an integer scalar
 *
 * @return R_NilValue
 */
SEXP mer_coefGets(SEXP x, SEXP coef, SEXP pType)
{
    int clen = coef_length(LENGTH(GET_SLOT(x, lme4_flistSym)),
			   INTEGER(GET_SLOT(x, lme4_ncSym)));
    if (LENGTH(coef) != clen || !isReal(coef))
	error(_("coef must be a numeric vector of length %d"), clen);
    internal_mer_coefGets(x, REAL(coef), asInteger(pType));
    return x;
}

/**
 * The naive way of calculating the trace of the hat matrix
 *
 * @param x pointer to an mer object
 *
 * @return trace of the hat matrix
 */

SEXP mer_hat_trace(SEXP x)
{
    SEXP Zt = GET_SLOT(x, lme4_ZtSym);
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    int *Zti = INTEGER(GET_SLOT(Zt, lme4_iSym)),
	*Ztp = INTEGER(GET_SLOT(Zt, lme4_pSym)), i, ione = 1, j,
	n = INTEGER(GET_SLOT(Zt, lme4_DimSym))[1],
	p = LENGTH(GET_SLOT(x, lme4_rXySym)),
	q = LENGTH(GET_SLOT(x, lme4_rZySym));
    double *Xcp = Calloc(n * p, double),
	*RXX = REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym), lme4_xSym)),
	*RZX = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXSym), lme4_xSym)),
	*Ztx = REAL(GET_SLOT(Zt, lme4_xSym)),
	*wrk = Calloc(q, double), m1 = -1, one = 1, tr;
    cholmod_dense *zrow = M_numeric_as_chm_dense(wrk, q);

    mer_factor(x);
    Memcpy(Xcp, REAL(GET_SLOT(x, lme4_XSym)), n * p);

    /* Accumulate F-norm of L^{-1}Zt and downdate rows of Xcp */
/* FIXME: Does this handle a non-trivial permutation properly? */
    for (j = 0, tr = 0; j < n; j++) { /* j'th column of Zt */
	cholmod_dense *sol; double *sx;
	for (i = 0; i < q; i++) wrk[i] = 0;
	for (i = Ztp[j]; i < Ztp[j + 1]; i++) wrk[Zti[i]] = Ztx[i];
	sol = M_cholmod_solve(CHOLMOD_L, L, zrow, &c);
	sx = (double*)(sol->x);
	for (i = 0; i < q; i++) tr += sx[i] * sx[i];
				/* downdate jth row of Xcp */
 	F77_CALL(dgemv)("T", &q, &p, &m1, RZX, &q, sx, &ione,
 			&one, Xcp + j, &n);
	M_cholmod_free_dense(&sol, &c);
    }
    F77_CALL(dtrsm)("R", "U", "N", "N", &n, &p, &one, RXX, &p, Xcp, &n);
    for (i = 0; i < n * p; i++) tr += Xcp[i] * Xcp[i];

    Free(zrow); Free(Xcp);
    return ScalarReal(tr);
}

/**
 * Faster calculation of the trace of the hat matrix (due to Jialiang Li)
 *
 * @param x pointer to an mer object
 *
 * @return trace of the hat matrix
 */

SEXP mer_hat_trace2(SEXP x)
{
    SEXP Omega = GET_SLOT(x, lme4_OmegaSym),
	ncp = GET_SLOT(x, lme4_ncSym);
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*nc = INTEGER(ncp),
	nf = LENGTH(ncp), i, j, k,
	p = LENGTH(GET_SLOT(x, lme4_rXySym)),
	q = LENGTH(GET_SLOT(x, lme4_rZySym));
    double
	*RZXicp = Calloc(q * p, double),
	one = 1, tr = p + q;
				/* factor and evaluate RZXinv */
    mer_factor(x);
    internal_mer_RZXinv(x);
    Memcpy(RZXicp, REAL(GET_SLOT(GET_SLOT(x, lme4_RZXinvSym),
				 lme4_xSym)), q * p);
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	int ncisqr = nci * nci, nlev = (Gp[i + 1] - Gp[i])/nci;
	double *deli = REAL(GET_SLOT(M_dpoMatrix_chol(VECTOR_ELT(Omega, i)),
				     lme4_xSym));
	cholmod_sparse *sol, *Prhs,
	    *rhs = M_cholmod_allocate_sparse(q, nci, ncisqr, TRUE, TRUE,
					     0, CHOLMOD_REAL, &c);
	int *rhi = (int *)(rhs->i), *rhp = (int *)(rhs->p);
	double *rhx = (double *)(rhs->x);

	rhp[0] = 0;		/* Establish column pointers and value of rhs */
	for (j = 0; j < nci; j++) {
	    rhp[j+1] = rhp[j] + nci;
	    for (k = 0; k < nci; k++) { /* transpose of deli */
		rhx[j * nci + k] = deli[k * nci + j];
		rhi[j * nci + k] = Gp[i] + k; /* initial row numbers */
	    }
	    /* zero the upper triangle (just in case) */
	    for (k = 0; k < j; k++) rhx[j * nci + k] = 0;
	}
	for (j = 0; j < nlev; j++) {
	    /* Evaluate nci rows of Delta RZXinv */
	    F77_CALL(dtrmm)("L", "U", "N", "N", &nci, &p, &one, deli, &nci,
			    RZXicp + Gp[i] + j * nci, &q);
	    /* Solve for and accumulate nci columns of L^{-1} P Delta' */
	    Prhs = M_cholmod_spsolve(CHOLMOD_P, L, rhs, &c);
	    sol = M_cholmod_spsolve(CHOLMOD_L, L, Prhs, &c);
	    M_cholmod_free_sparse(&Prhs, &c);
	    for (k = 0; k < ((int*)(sol->p))[nci]; k++) {
		double sxk = ((double*)(sol->x))[k];
		tr -= sxk * sxk;
	    }
	    M_cholmod_free_sparse(&sol, &c);
	    /* Update rhs for the next set of rows */
	    for (k = 0; k < ncisqr; k++) rhi[k] += nci;
	}
	M_cholmod_free_sparse(&rhs, &c);
    }
    for (i = 0; i < q * p; i++) tr -= RZXicp[i] * RZXicp[i];
    Free(RZXicp);
    return ScalarReal(tr);
}

/**
 * Return L as a dtCMatrix object
 *
 * @param x pointer to an mer object
 *
 * @return L as an dtCMatrix object
 */
SEXP mer_dtCMatrix(SEXP x)
{
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym)), *Lcp;
    cholmod_sparse *Lm;
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    int *dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2)),
	nz, q;

    dims[0] = dims[1] = q = (int)(L->n);
    Lcp = M_cholmod_copy_factor(L, &c); Free(L); /* next call changes Lcp */
    Lm = M_cholmod_factor_to_sparse(Lcp, &c); M_cholmod_free_factor(&Lcp, &c);
    SET_SLOT(ans, lme4_uploSym, mkString("L"));
    SET_SLOT(ans, lme4_diagSym, mkString("N"));
    Memcpy(INTEGER(ALLOC_SLOT(ans, lme4_pSym, INTSXP, q + 1)),
	   (int*) Lm->p, q + 1);
    nz = ((int*)(Lm->p))[q];
    Memcpy(INTEGER(ALLOC_SLOT(ans, lme4_iSym, INTSXP, nz)),
	   (int*) Lm->i, nz);
    Memcpy(REAL(ALLOC_SLOT(ans, lme4_xSym, REALSXP, nz)),
	   (double*) Lm->x, nz);
    M_cholmod_free_sparse(&Lm, &c);
    UNPROTECT(1);
    return ans;
}

/**
 * Return L^{-1} as a dtCMatrix object
 *
 * @param x pointer to an mer object
 *
 * @return L^{-1} as an dtCMatrix object
 */
SEXP mer_dtCMatrix_inv(SEXP x)
{
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    cholmod_sparse
	*b = M_cholmod_allocate_sparse(L->n, L->n, L->n, 1, 1,
				       0, CHOLMOD_REAL, &c),
	*Linv;
    double *bx = (double *)(b->x);
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    int *bi = (int *) (b->i), *bp = (int *) (b->p),
	*dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2)),
	j, nz, q;

    dims[0] = dims[1] = q = (int)(L->n);
    for (j = 0; j < q; j++) {
	bp[j] = bi[j] = j;
	bx[j] = 1;
    }
    bp[q] = q;
    Linv = M_cholmod_spsolve(CHOLMOD_L, L, b, &c);
    M_cholmod_free_sparse(&b, &c);
    SET_SLOT(ans, lme4_uploSym, mkString("L"));
    SET_SLOT(ans, lme4_diagSym, mkString("N"));
    Memcpy(INTEGER(ALLOC_SLOT(ans, lme4_pSym, INTSXP, q + 1)),
	   (int *) Linv->p, q + 1);
    nz = ((int *)(Linv->p))[q];
    Memcpy(INTEGER(ALLOC_SLOT(ans, lme4_iSym, INTSXP, nz)),
	   (int *) Linv->i, nz);
    Memcpy(REAL(ALLOC_SLOT(ans, lme4_xSym, REALSXP, nz)),
	   (double *) Linv->x, nz);
    M_cholmod_free_sparse(&Linv, &c);
    UNPROTECT(1);
    return ans;
}

/**
 * Create and factor Z'Z+Omega if it has not already.
 * Also create RZX and RXX, the deviance components,
 * and the value of the deviance for both ML and  REML.
 *
 * @param x pointer to an lmer object
 *
 * @return NULL
 */
SEXP mer_factor(SEXP x)
{
    int *status = INTEGER(GET_SLOT(x, lme4_statusSym));
    if (!status[0]) {
	SEXP rXyP = GET_SLOT(x, lme4_rXySym),
	    rZyP = GET_SLOT(x, lme4_rZySym);
	int i, info, ione = 1, p = LENGTH(rXyP), q = LENGTH(rZyP);
	cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
	double *RXX = REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym), lme4_xSym)),
	    *RZX = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXSym), lme4_xSym)),
	    *rXy = REAL(rXyP), *rZy = REAL(rZyP),
	    *dcmp = REAL(GET_SLOT(x, lme4_devCompSym)),
	    *dev = REAL(GET_SLOT(x, lme4_devianceSym)),
	    one[2] = {1, 0}, m1[2] = {-1, 0};
	double nml = dcmp[0], nreml = dcmp[0] - dcmp[1];

	/* Inflate Z'Z to Z'Z+Omega and factor to form L. Form RZX and
	 * rZy. Update stage flag, dcmp[4] and dcmp[5]. */
	internal_mer_Zfactor(x, L);
				/* downdate XtX and factor */
	if ((info = internal_mer_Xfactor(x))) /* unable to factor downdated XtX */
	    error(_("Leading minor of order %d in downdated X'X is not positive definite"),
		  info);
	for (dcmp[6] = 0, i = 0; i < p; i++) /* 2 * logDet(RXX) */
	    dcmp[6] += 2. * log(RXX[i * (p + 1)]);
				/* solve for rXy  and ryy^2 */
	Memcpy(rXy, REAL(GET_SLOT(x, lme4_XtySym)), p);
	F77_CALL(dgemv)("T", &q, &p, m1, RZX, &q, rZy, &ione, one, rXy, &ione);
	F77_CALL(dtrsv)("U", "T", "N", &p, RXX, &p, rXy, &ione);
	dcmp[3] = log(dcmp[2] /* dcmp[3] = log(ryy^2); dcmp[2] = y'y; */
		      - F77_CALL(ddot)(&p, rXy, &ione, rXy, &ione)
		      - F77_CALL(ddot)(&q, rZy, &ione, rZy, &ione));
				/* evaluate ML and REML deviance */
	dev[0] = dcmp[4] - dcmp[5] +
	    nml*(1.+dcmp[3]+log(2.*PI/nml));
	dev[1] = dcmp[4] - dcmp[5] + dcmp[6] +
	    nreml*(1.+dcmp[3]+log(2.*PI/nreml));
	if (dcmp[7] >= 0) dcmp[7] = internal_mer_sigma(x, -1);
	Free(L);
    }
    return R_NilValue;
}

/**
 * Return the fitted values as an SEXP
 *
 * @param x pointer to an mer object
 * @param useFe pointer to a logical scalar indicating if the fixed
 * effects should be used
 * @param useRe pointer to a logical scalar indicating if the random
 * effects should be used
 *
 * @return pointer to a numeric array of fitted values
 */

SEXP mer_fitted(SEXP x)
{
    int n = LENGTH(GET_SLOT(x, lme4_ySym));
    SEXP ans = PROTECT(allocVector(REALSXP, n));

    internal_mer_fitted(x, (double*) NULL, REAL(ans));
    UNPROTECT(1);
    return ans;
}

/**
 * Extract the conditional estimates of the fixed effects
 *
 * @param x Pointer to an mer object
 *
 * @return a numeric vector containing the conditional estimates of
 * the fixed effects
 */
SEXP mer_fixef(SEXP x)
{
    int nf = LENGTH(GET_SLOT(x, lme4_OmegaSym));
    SEXP ans;

    mer_secondary(x);
    ans = PROTECT(duplicate(GET_SLOT(x, lme4_fixefSym)));
    setAttrib(ans, R_NamesSymbol,
	      duplicate(VECTOR_ELT(GET_SLOT(x, lme4_cnamesSym), nf)));
    UNPROTECT(1);
    return ans;
}

/**
 * Fill in the gradComp and bVar slots.  Each component in the gradComp slot
 * consists of four symmetric matrices used to generate the gradient or the ECME
 * step.  They are
 *  1) -m_i\bOmega_i^{-1}
 *  2) \bB_i\bB_i\trans
 *  3) \tr\left[\der_{\bOmega_i}\bOmega\left(\bZ\trans\bZ+\bOmega\right)\inv\right]
 *  4) The term added to 3) to get \tr\left[\der_{\bOmega_i}\bOmega\vb\right]
 *
 * @param x pointer to an mer object
 * @param val pointer to a list of matrices of the correct sizes
 *
 * @return NULL
 */
SEXP mer_gradComp(SEXP x)
{
    int *status = INTEGER(GET_SLOT(x, lme4_statusSym));

    mer_secondary(x);
    if (status[0] < 3) {
	SEXP bVarP = GET_SLOT(x, lme4_bVarSym),
	    OmegaP = GET_SLOT(x, lme4_OmegaSym),
	    gradComp = GET_SLOT(x, lme4_gradCompSym),
	    ranefP = GET_SLOT(x, lme4_ranefSym);
	int q = LENGTH(ranefP), p = LENGTH(GET_SLOT(x, lme4_rXySym));
	cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
	int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	    *nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	    i, j, k, nf = length(OmegaP);
	double *b = REAL(GET_SLOT(x, lme4_ranefSym)),
	    *RZXinv = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXinvSym),
				    lme4_xSym)),
	    alpha;

	alpha = 1./internal_mer_sigma(x, -1);
	alpha = alpha * alpha;

	internal_mer_RZXinv(x);
	internal_mer_bVar(x);
	for (i = 0; i < nf; i++) {
	    int nci = nc[i], RZXrows = Gp[i + 1] - Gp[i];
	    int ncisq = nci * nci, nlev = RZXrows/nci;
	    double *bVi = REAL(VECTOR_ELT(bVarP, i)),
		*bi = b + Gp[i], *mm = REAL(VECTOR_ELT(gradComp, i)),
		*tmp = Memcpy(Calloc(ncisq, double),
			      REAL(GET_SLOT(M_dpoMatrix_chol(VECTOR_ELT(OmegaP, i)),
					    lme4_xSym)), ncisq),
		*RZXi = RZXinv + Gp[i],
		dlev = (double) nlev,
		one[] = {1,0}, zero[] = {0,0};

	    if (nci == 1) {
		int ione = 1;
		mm[0] = ((double) nlev)/(tmp[0] * tmp[0]);
		mm[1] = alpha * F77_CALL(ddot)(&nlev, bi, &ione, bi, &ione);
		mm[2] = 0.;
		for (k = 0; k < nlev; k++) mm[2] += bVi[k];
		mm[3] = 0.;
		for (j = 0; j < p; j++) {
		    mm[3] += F77_CALL(ddot)(&RZXrows, RZXi + j * q, &ione,
					    RZXi + j * q, &ione);
		}
	    } else {
		AZERO(mm, 4 * ncisq);
		F77_CALL(dtrtri)("U", "N", &nci, tmp, &nci, &j);
		if (j)
		    error(_("Omega[[%d]] is not positive definite"), i + 1);
		F77_CALL(dsyrk)("U", "N", &nci, &nci, &dlev, tmp, &nci,
				zero, mm, &nci);
		mm += ncisq;	/* \bB_i term */
		F77_CALL(dsyrk)("U", "N", &nci, &nlev, &alpha, bi, &nci,
				zero, mm, &nci);
		mm += ncisq;     /* Sum of diagonal blocks of the inverse
				  * (Z'Z+Omega)^{-1} */
		for (j = 0; j < ncisq; j++) {
		    for (k = 0; k < nlev; k++) mm[j] += bVi[j + k*ncisq];
		}
		mm += ncisq;	/* Extra term for \vb */
		for (j = 0; j < p; j++) {
		    F77_CALL(dsyrk)("U", "N", &nci, &nlev, one,
				    RZXi + j * q, &nci,
				    one, mm, &nci);
		}
	    }
	    Free(tmp);
	}
	Free(L);
	status[0] = 3;
    }
    return R_NilValue;
}

/**
 * Evaluate the gradient vector
 *
 * @param x Pointer to an lmer object
 * @param pType Pointer to an integer indicator of the
 * parameterization being used
 *
 * @return pointer to a gradient vector
 */
SEXP mer_gradient(SEXP x, SEXP pType)
{
    SEXP Omega = GET_SLOT(x, lme4_OmegaSym);
    SEXP gradComp = GET_SLOT(x, lme4_gradCompSym);
    int *nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	dind, i, ifour = 4, ione = 1, nf = length(Omega),
	odind, ptyp = asInteger(pType);
    int REML = INTEGER(GET_SLOT(x, lme4_statusSym))[1];
    SEXP val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));
    double cc[] = {-1, 1, 1, REML ? 1 : 0},
      	*valp = REAL(val),
	one = 1.0, zero = 0.0;

    mer_gradComp(x);
    dind = 0;			/* index into val for diagonals */
    for (i = 0; i < nf; i++) {
	SEXP Omgi = VECTOR_ELT(Omega, i);
	int nci = nc[i], ncisqr = nci * nci;
	double *tmp = Calloc(ncisqr, double);

	F77_CALL(dgemm)("N", "N", &ncisqr, &ione, &ifour, &one,
			REAL(VECTOR_ELT(gradComp, i)), &ncisqr,
			cc, &ifour, &zero, tmp, &ncisqr);
	if (nci == 1) {
	    double omg = *REAL(GET_SLOT(Omgi, lme4_xSym));
	    valp[dind++] =
		(ptyp?((ptyp == 1)? omg : -omg * omg) : 1) * tmp[0];
	} else {
	    int ii, j, ncip1 = nci + 1;

	    odind = dind + nci; /* index into val for off-diagonals */
	    if (ptyp) {
		double *chol = REAL(GET_SLOT(M_dpoMatrix_chol(Omgi), lme4_xSym)),
		    *tmp2 = Calloc(ncisqr, double);

		/* Overwrite the gradient with respect to positions in
		 * Omega[[i]] by the gradient with respect to the
		 * unconstrained parameters.*/

		/* tmp2 := chol %*% tmp using only upper triangle of tmp */
		F77_CALL(dsymm)("R", "U", &nci, &nci, &one, tmp, &nci,
				chol, &nci, &zero, tmp2, &nci);
		/* full symmetric product gives diagonals */
		F77_CALL(dtrmm)("R", "U", "T", "N", &nci, &nci, &one, chol,
				&nci, Memcpy(tmp, tmp2, ncisqr), &nci);
		/* overwrite upper triangle with gradients for L' */
		for (ii = 1; ii < nci; ii++) {
		    for (j = 0; j < ii; j++) {
			tmp[j + ii*nci] = chol[j*ncip1] * tmp2[j + ii*nci];
			tmp[ii + j*nci] = 0.;
		    }
		}
		if (ptyp > 1)
		    for (ii = 0; ii < nci; ii++) {
			int ind = ii * ncip1;
			double sqrtd = chol[ind];
			tmp[ind] *= -(sqrtd*sqrtd);
		    }
		Free(tmp2);
	    }
	    for (j = 0; j < nci; j++) {
		valp[dind + j] = tmp[j * ncip1];
		for (ii = 0; ii < j; ii++) /* offdiagonals count twice */
		    valp[odind++] = 2. * tmp[ii + j * nci];
	    }
	    dind = odind;
	}
	Free(tmp);
    }
    UNPROTECT(1);
    return val;
}

/**
 * Create and insert initial values for Omega.
 *
 * @param x pointer to an mer object
 *
 * @return NULL
 */
SEXP mer_initial(SEXP x)
{
    SEXP Omg = GET_SLOT(x, lme4_OmegaSym),
	ZtZ = GET_SLOT(x, lme4_ZtZSym);
    int	*Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	*p = INTEGER(GET_SLOT(ZtZ, lme4_pSym)),
	i, nf = length(Omg);
    double *xx = REAL(GET_SLOT(ZtZ, lme4_xSym));

    for (i = 0; i < nf; i++) {
	SEXP Omgi = VECTOR_ELT(Omg, i);
	double *omgi = REAL(GET_SLOT(Omgi, lme4_xSym));
	int bb = Gp[i], j, k, nci = nc[i];
	int ncip1 = nci + 1, nlev = (Gp[i + 1] - bb)/nci;

	AZERO(omgi, nci * nci);
	for (j = 0; j < nlev; j++) {
	    int base = bb + j * nci;
	    for (k = 0; k < nci; k++)
		/* add the last element in the column */
		omgi[k * ncip1] += xx[p[base + k + 1] - 1];
	}
	for (k = 0; k < nci; k++) omgi[k * ncip1] *= 0.375/nlev;
	SET_SLOT(Omgi, lme4_factorSym, allocVector(VECSXP, 0));
	M_dpoMatrix_chol(Omgi);
    }
    flag_not_factored(x);
    return R_NilValue;
}

/**
 * Externally callable check on nesting
 *
 * @param x Pointer to an mer object
 *
 * @return a scalar logical value indicating if ZtZ corresponds to a
 * simple nested structure.
 */
SEXP mer_isNested(SEXP x)
{
    cholmod_sparse *ZtZ = M_as_cholmod_sparse(GET_SLOT(x, lme4_ZtZSym));
    cholmod_sparse *ZtZl = M_cholmod_transpose(ZtZ, (int) ZtZ->xtype, &c);
    SEXP ncp = GET_SLOT(x, lme4_ncSym);
    int ans = internal_mer_isNested(LENGTH(ncp), INTEGER(ncp),
				    INTEGER(GET_SLOT(x, lme4_GpSym)),
				    ZtZl->p);
    Free(ZtZ); M_cholmod_free_sparse(&ZtZl, &c);
    return ScalarLogical(ans);
}

/**
 * Extract the conditional modes of the random effects.
 *
 * @param x Pointer to an mer object
 *
 * @return a list of matrices containing the conditional modes of the
 * random effects
 */
SEXP mer_ranef(SEXP x)
{
    SEXP cnames = GET_SLOT(x, lme4_cnamesSym),
	flist = GET_SLOT(x, lme4_flistSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	i, ii, jj,
	nf = length(flist);
    SEXP val = PROTECT(allocVector(VECSXP, nf));
    double *b = REAL(GET_SLOT(x, lme4_ranefSym));

    mer_secondary(x);
    setAttrib(val, R_NamesSymbol,
	      duplicate(getAttrib(flist, R_NamesSymbol)));
    for (i = 0; i < nf; i++) {
	SEXP nms, rnms = getAttrib(VECTOR_ELT(flist, i), R_LevelsSymbol);
	int nci = nc[i], mi = length(rnms);
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
 * Update the secondary slots - fixef and ranef
 *
 * @param x pointer to a mer object
 *
 */
SEXP mer_secondary(SEXP x)
{
    int *status = INTEGER(GET_SLOT(x, lme4_statusSym));

    mer_factor(x);
    if (status[0] < 2) {
	internal_mer_fixef(x);
	internal_mer_ranef(x);
    }
    return R_NilValue;
}

/**
 * Extract the posterior variances of the random effects
 *
 * @param x pointer to a mer object
 *
 * @return pointer to a list of arrays
 */
SEXP mer_postVar(SEXP x)
{
    SEXP ans;
    double *RZXi = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXinvSym), lme4_xSym)),
	*dcmp = REAL(GET_SLOT(x, lme4_devCompSym)), one = 1, sc;
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	i, ione = 1, nf,
	p = LENGTH(GET_SLOT(x, lme4_rXySym)),
	q = LENGTH(GET_SLOT(x, lme4_rZySym));

    sc = dcmp[7] * dcmp[7];
    mer_gradComp(x);
    ans = PROTECT(duplicate(GET_SLOT(x, lme4_bVarSym)));
    nf = LENGTH(ans);
    for (i = 0; i < nf; i++) {
	SEXP ansi = VECTOR_ELT(ans, i);
	int *dims = INTEGER(getAttrib(ansi, R_DimSymbol));
	int j, nc = dims[1], nlev = dims[2];
	int ncsqr = dims[0] * nc;
	int ntot = ncsqr * nlev;
	double *vv = REAL(ansi);

	if (dims[0] != nc)
	    error(_("rows and columns of element %d of bVar do not match"),
		  i + 1);
	for (j = 0; j < nlev; j++)
	    F77_CALL(dsyrk)("U", "N", &nc, &p,
			    &one, RZXi + Gp[i] + j * nc, &q,
			    &one, vv + j * ncsqr, &nc);
	if (sc != 1) F77_CALL(dscal)(&ntot, &sc, vv, &ione);
	if (nc > 1) {
	    for (j = 0; j < nlev; j++)
		internal_symmetrize(vv + j * ncsqr, nc);
	}
    }
    UNPROTECT(1);
    return ans;
}

/**
 * Extract the ML or REML conditional estimate of sigma
 *
 * @param x pointer to an mer object
 * @param REML logical scalar - TRUE if REML estimates are requested
 *
 * @return pointer to a numeric scalar
 */
SEXP mer_sigma(SEXP x, SEXP REML)
{
    return ScalarReal(
	internal_mer_sigma(x,
			   (REML == R_NilValue) ? -1 :
			   (asLogical(REML))));
}

/**
 * Simulate a set of linear predictors from the random effects part of
 * an mer object
 *
 * @param x Pointer to an mer object
 * @param np Pointer to an integer giving the number of values to simulate
 *
 * @return a matrix of simulated linear predictors
 */
SEXP mer_simulate(SEXP x, SEXP nsimP)
{
    int *nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	*Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	i, ii, j, nsim = asInteger(nsimP),
	nf = LENGTH(GET_SLOT(x, lme4_OmegaSym)),
	n = LENGTH(GET_SLOT(x, lme4_ySym)),
	q = LENGTH(GET_SLOT(x, lme4_ZtySym));
    SEXP ans = PROTECT(allocMatrix(REALSXP, n, nsim)),
	Omega = GET_SLOT(x, lme4_OmegaSym);
    cholmod_dense *cha = M_as_cholmod_dense(ans),
	*chb = M_cholmod_allocate_dense(q, nsim, q, CHOLMOD_REAL, &c);
    double *dcmp = REAL(GET_SLOT(x, lme4_devCompSym)),
	one[] = {1,0}, zero[] = {0,0};
    double scale = (dcmp[7] < 0) ? -dcmp[7] : dcmp[7];
    cholmod_sparse *Zt = M_as_cholmod_sparse(GET_SLOT(x, lme4_ZtSym));

    GetRNGstate();
    for (ii = 0; ii < nsim; ii++) {
	for (i = 0; i < nf; i++) {
	    int nci = nc[i], relen = Gp[i + 1] - Gp[i];
	    int nlev = relen/nci;
	    double *bi = (double *)(chb->x) + ii * q + Gp[i],
		*Rmat = REAL(GET_SLOT(M_dpoMatrix_chol(VECTOR_ELT(Omega, i)),
				      lme4_xSym));

	    for (j = 0; j < relen; j++) bi[j] = norm_rand();
	    F77_CALL(dtrsm)("L", "U", "N", "N", &nci, &nlev, &scale,
			    Rmat, &nci, bi, &nci);
	}
    }
    PutRNGstate();

    if (!M_cholmod_sdmult(Zt, 1, one, zero, chb, cha, &c))
	error(_("cholmod_sdmult failed"));
    M_cholmod_free_dense(&chb, &c);
    Free(Zt); Free(cha);
    UNPROTECT(1);
    return ans;
}

/**
 * Externally callable version of internal_mer_update_ZXy
 *
 * @param x pointer to an mer object
 *
 * @return NULL
 */
SEXP mer_update_ZXy(SEXP x)
{
    internal_mer_update_ZXy(x,
			    INTEGER(GET_SLOT(GET_SLOT(x, lme4_LSym),
					     lme4_permSym)));
    return R_NilValue;
}

/**
 * Update the y slot (and slots derived from it) in an mer object
 *
 * @param x pointer to an mer object
 * @param ynew pointer to a numeric vector of length n
 *
 * @return NULL
 */
SEXP mer_update_y(SEXP x, SEXP ynew)
{
    SEXP y = GET_SLOT(x, lme4_ySym),
	Xty = GET_SLOT(x, lme4_XtySym),
	Zty = GET_SLOT(x, lme4_ZtySym);
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    int *perm = (int*)(L->Perm), i, ione = 1,
	n = LENGTH(y), p = LENGTH(Xty), q = LENGTH(Zty);
    cholmod_sparse *Zt = M_as_cholmod_sparse(GET_SLOT(x, lme4_ZtSym));
    cholmod_dense *td1, *yd = M_as_cholmod_dense(GET_SLOT(x, lme4_ySym));
    double one[] = {1,0}, zero[] = {0,0};

    if (!isReal(ynew) || LENGTH(ynew) != n)
	error(_("ynew must be a numeric vector of length %d"), n);
    Memcpy(REAL(y), REAL(ynew), n);
    				/* y'y */
    REAL(GET_SLOT(x, lme4_devCompSym))[2] =
	F77_CALL(ddot)(&n, REAL(y), &ione, REAL(y), &ione);
				/* PZ'y into Zty */
    td1 = M_cholmod_allocate_dense(q, 1, q, CHOLMOD_REAL, &c);
    if (!M_cholmod_sdmult(Zt, 0, one, zero, yd, td1, &c))
	error(_("cholmod_sdmult failed"));
    for (i = 0; i < q; i++) REAL(Zty)[i] = ((double *)(td1->x))[perm[i]];
    M_cholmod_free_dense(&td1, &c); Free(yd); Free(Zt);
    				/* Xty */
    F77_CALL(dgemv)("T", &n, &p, one, REAL(GET_SLOT(x, lme4_XSym)),
		    &n, REAL(y), &ione, zero, REAL(Xty), &ione);
    flag_not_factored(x);
    Free(L);
    return R_NilValue;
}

/**
 * Check validity of an mer object.
 *
 * @param x Pointer to an mer object
 *
 * @return TRUE if the object is a valid lmer object, else a string
 * describing the nature of the violation.
 */
SEXP mer_validate(SEXP x)
{
    SEXP
	ZtXP = GET_SLOT(x, lme4_ZtXSym),
	XtXP = GET_SLOT(x, lme4_XtXSym),
	RZXP = GET_SLOT(x, lme4_RZXSym),
	RXXP = GET_SLOT(x, lme4_RXXSym)
	/* , cnames = GET_SLOT(x, lme4_cnamesSym) */
	;
    int *ZtXd = INTEGER(getAttrib(ZtXP, lme4_DimSym)),
	*XtXd = INTEGER(getAttrib(XtXP, lme4_DimSym));

    if (!match_mat_dims(ZtXd, INTEGER(getAttrib(RZXP, lme4_DimSym))))
	return mkString(_("Dimensions of slots ZtX and RZX must match"));
    if (!match_mat_dims(XtXd, INTEGER(getAttrib(RXXP, lme4_DimSym))))
	return mkString(_("Dimensions of slots XtX and RXX must match"));
    if (ZtXd[1] != XtXd[0] || XtXd[0] != XtXd[1])
	return mkString(_("Slot XtX must be a square matrix with same ncol as ZtX"));
    return ScalarLogical(1);
}

static
SEXP SEXP_Zt(int n, int ii, SEXP fi, SEXP tmmat)
{
    int *dims = INTEGER(getAttrib(tmmat, R_DimSymbol)), *fac =INTEGER(fi), j, k;
    int m = dims[0], nlev = LENGTH(getAttrib(fi, R_LevelsSymbol));
    SEXP ans = PROTECT(alloc_dgCMatrix(m * nlev, n, m * n, R_NilValue, R_NilValue));
    int *i = INTEGER(GET_SLOT(ans, lme4_iSym)), *p = INTEGER(GET_SLOT(ans, lme4_pSym));
    
    if (!isFactor(fi) || LENGTH(fi) != n)
	error(_("fl[[%d]] must be a factor of length %d"), ii + 1, n);
    if (!isMatrix(tmmat) || !isReal(tmmat))
	error(_("Ztl[[%d]] must be real matrix"), ii + 1);
    if (dims[1] != n)
	error(_("Ztl[[%d]] must have %d columns"), ii + 1, n);

    p[n] = m * n;
    for (j = 0; j < n; j++) {
	p[j] = m * j;
	for (k = 0; k < m; k++) i[j * m + k] = (fac[j] - 1) * m + k;
    }
    Memcpy(REAL(GET_SLOT(ans, lme4_xSym)), REAL(tmmat), m * n);
    UNPROTECT(1);
    return ans;
}

/**
 * Create a list of sparse Zt matrices from a factor list and a list
 * of dense, skinny model matrices
 *
 * @param fl list of factors
 * @param Ztl list of transposes of model matrices
 *
 * @return a list of sparse (full) Zt matrices
 */
SEXP Ztl_sparse(SEXP fl, SEXP Ztl)
{
    int i, nf = LENGTH(fl), nobs = LENGTH(VECTOR_ELT(fl, 0));
    SEXP ans = PROTECT(allocVector(VECSXP, nf));
    
    setAttrib(ans, R_NamesSymbol, duplicate(getAttrib(fl, R_NamesSymbol)));
    for (i = 0; i < nf; i++)
	SET_VECTOR_ELT(ans, i, SEXP_Zt(nobs, i, VECTOR_ELT(fl, i), VECTOR_ELT(Ztl, i)));
    UNPROTECT(1);
    return ans;
}

/**
 * Create a new sparse Zt matrix by carrying over elements for the same level of f
 *
 * @param f factor determining the carryover (e.g. student)
 * @param Zt sparse model matrix for another factor (e.g. teacher)
 *
 * @return modified model matrix
 */
SEXP Zt_carryOver(SEXP fp, SEXP Zt)
{
    cholmod_sparse *ans, *chsz = M_as_cholmod_sparse(Zt);
    cholmod_triplet *ant, *chtz = M_cholmod_sparse_to_triplet(chsz, &c);
    int *cct, *p = (int*)(chsz->p), *f = INTEGER(fp);
    int cmax, j, jj, k, last, n = LENGTH(fp), nlev, nnz, ntot, q = p[1] - p[0];
    int *ii, *ij, *oi, *oj, ip, op;
    double *ix, *ox;

    Free(chsz);
    if (!isFactor(fp)) error(_("f must be a factor"));
    nlev = LENGTH(getAttrib(fp, R_LevelsSymbol));
    cct = Calloc(nlev, int);
    
    if (chtz->ncol != n) error(_("ncol(Zt) must match length(fp)"));
    for (j = 0; j < n; j++)	/* check consistency of p */
	if (p[j+1] - p[j] != q) error(_("nonzeros per column in Zt must be constant"));
				/* create column counts */
    for (last = -1, j = 0; j < n; j++) {
	int ll = f[j] - 1;
	if (last > ll) error(_("f is not in increasing order"));
	if (ll == last) cct[ll]++;
	else {cct[ll] = 1; last = ll;};
    }
    if (nlev != f[n - 1]) error(_("number of levels of f is not consistent"));
				/*  max column count and total nonzeros*/
    for (cmax = -1, ntot = 0, k = 0; k < nlev; k++) {
	cmax = (cct[k] > cmax) ? cct[k] : cmax;
	ntot += (cct[k] * (cct[k] + 1))/2;
    }
    nnz = ntot * q;
    ant = M_cholmod_allocate_triplet(chtz->nrow, chtz->ncol, (size_t)(nnz), 
				     0 /*stype*/, CHOLMOD_REAL, &c);
    ip = 0; ii = (int*)(chtz->i); ij = (int*)(chtz->j); ix = (double*)(chtz->x);
    op = 0; oi = (int*)(ant->i); oj = (int*)(ant->j); ox = (double*)(ant->x);
    for (k = 0; k < nlev; k++) {
	for (j = 0; j < cct[k]; j++) {
	    for (jj = 0; jj < cct[k] - j; jj++) {
		oi[op] = ii[ip]; oj[op] = ij[ip] + jj; ox[op] = ix[ip];
		op++;
	    }
	    ip++;
	}
    }
    ant->nnz = nnz;
#if 0
				/* fake it for the time being */
    nnz = ant->nnz = chtz->nnz;
    Memcpy((int*)(ant->i), (int*)(chtz->i), nnz);
    Memcpy((int*)(ant->j), (int*)(chtz->j), nnz);
    Memcpy((double*)(ant->x), (double*)(chtz->x), nnz);
#endif
    ans = M_cholmod_triplet_to_sparse(ant, nnz, &c);
    M_cholmod_free_triplet(&chtz, &c);
    M_cholmod_free_triplet(&ant, &c);
    Free(cct);
    return M_chm_sparse_to_SEXP(ans, 1, 0, 0, "", R_NilValue);
}
	
#if 0
static
cholmod_sparse *chm_Zt(int n, int ii, SEXP fi, SEXP tmmat)
{
    cholmod_sparse *ans;
    int *dims, *fac, *i, *p, j, k, m;

    if (!isFactor(fi) || LENGTH(fi) != n)
	error(_("fl[[%d]] must be a factor of length %d"), ii + 1, n);
    if (!isMatrix(tmmat) || !isReal(tmmat))
	error(_("Ztl[[%d]] must be real matrix"), ii + 1);
    dims = INTEGER(getAttrib(tmmat, R_DimSymbol));
    if (dims[1] != n)
	error(_("Ztl[[%d]] must have %d columns"), ii + 1, n);
    m = dims[0];
    ans = M_cholmod_allocate_sparse(m * LENGTH(getAttrib(fi, R_LevelsSymbol)),
				    n, m * n, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
    i = (int *)(ans->i); p = (int *)(ans->p); fac = INTEGER(fi);

    p[n] = m * n;
    for (j = 0; j < n; j++) {
	p[j] = m * j;
	for (k = 0; k < m; k++) i[j * m + k] = (fac[j] - 1) * m + k;
    }
    Memcpy((double*)(ans->x), REAL(tmmat), m * n);
    return ans;
}

/**
 * Create the sparse Zt matrix from a factor list and list of model matrices
 *
 * @param fl list of factors
 * @param Ztl list of transposes of model matrices
 *
 * @return a freshly created sparse Zt object
 */
SEXP Zt_create1(SEXP fl, SEXP Ztl)
{
    cholmod_sparse *ans, *cat, *cur;
    SEXP f0;
    int i, nf = LENGTH(fl), nobs;

    if (!isNewList(fl) || nf < 1)
	error(_("fl must be a non-null list"));
    if (!isNewList(Ztl) || LENGTH(Ztl) != nf)
	error(_("Ztl must be a list of length %d"), nf);
    f0 = VECTOR_ELT(fl, 0);
    nobs = LENGTH(f0);
    if (!isFactor(f0) || nobs < 1)
	error(_("fl[[1]] must be a non-empty factor"));
    ans = chm_Zt(nobs, 0, f0, VECTOR_ELT(Ztl, 0));
    for (i = 1; i < nf; i++) {	/* check consistency, create and rbind a layer */
	cur = chm_Zt(nobs, i, VECTOR_ELT(fl, i), VECTOR_ELT(Ztl, i));
	cat = M_cholmod_vertcat(ans, cur, TRUE, &c);
	M_cholmod_free_sparse(&ans, &c); M_cholmod_free_sparse(&cur, &c);
	ans = cat;
    }
    return M_chm_sparse_to_SEXP(ans, 1, 0, 0, "", R_NilValue);
}

/**
 * Create the sparse Zt matrix with carry-over from a factor list and
 * list of model matrices
 *
 * @param fl list of factors
 * @param Ztl list of transposes of model matrices
 *
 * @return a freshly created sparse Zt object
 */
SEXP Zt_create(SEXP fl, SEXP Ztl)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix"))), fi, tmmat;
    int *dims, *p, *ii, i, nrtot = 0, nf = LENGTH(fl), nobs;
    int *Gp = Calloc(nf + 1, int), *nr = Calloc(nf, int),
	*offset = Calloc(nf + 1, int);
    double *x;

    if (!isNewList(fl) || nf < 1)
	error(_("fl must be a non-null list"));
    if (!isNewList(Ztl) || LENGTH(Ztl) != nf)
	error(_("Ztl must be a list of length %d"), nf);
    fi = VECTOR_ELT(fl, 0);
    nobs = LENGTH(fi);
    if (!isFactor(fi) || nobs < 1)
	error(_("fl[[1]] must be a non-empty factor"));
    offset[0] = Gp[0] = 0;
    for (i = 0; i < nf; i++) {	/* check consistency and establish dimensions */
	fi = VECTOR_ELT(fl, i);	/* grouping factor */
	if (!isFactor(fi) || LENGTH(fi) != nobs)
	    error(_("fl[[%d]] must be a factor of length %d"), i + 1, nobs);
	tmmat = VECTOR_ELT(Ztl, i); /* transpose of model matrix */
	if (!isMatrix(tmmat) || !isReal(tmmat))
	    error(_("Ztl[[%d]] must be real matrix"), i + 1);
	dims = INTEGER(getAttrib(tmmat, R_DimSymbol));
	if (dims[1] != nobs)
	    error(_("Ztl[[%d]] must have %d columns"), i + 1, nobs);
	nrtot += (nr[i] = dims[0]);
	offset[i + 1] = offset[i] + nr[i];
	Gp[i + 1] = Gp[i] + dims[0] * LENGTH(getAttrib(fi, R_LevelsSymbol));
    }
    dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2));
    dims[0] = Gp[nf]; dims[1] = nobs;
    p = INTEGER(ALLOC_SLOT(ans, lme4_pSym, INTSXP, nobs + 1));
    ii = INTEGER(ALLOC_SLOT(ans, lme4_iSym, INTSXP, nrtot * nobs));
    x = REAL(ALLOC_SLOT(ans, lme4_xSym, REALSXP, nrtot * nobs));
    p[0] = 0; for(i = 0; i < nobs; i++) p[i + 1] = p[i] + nrtot;

    for (i = 0; i < nf; i++) {	/* fill ans */
	int *vals = INTEGER(VECTOR_ELT(fl, i)), j;
	double *xvals = REAL(VECTOR_ELT(Ztl, i));

	for (j = 0; j < nobs; j++) {
	    int jbase = Gp[i] + nr[i] * (vals[j] - 1), k;
	    for (k = 0; k < nr[i]; k++) {
		int ind = j * nrtot + offset[i] + k;
		ii[ind] = jbase + k;
		x[ind] = xvals[j * nr[i] + k];
	    }
	}
    }

    Free(offset); Free(Gp); Free(nr);
    UNPROTECT(1);
    return ans;
}
#endif

/* MCMCsamp for generalized linear mixed models.
 *
 * 1) Detect if there  is a scale factor or not (not done yet).
 * 2) Sample beta and b according to the normal approximation
 * 3) Evaluate the Metropolis-Hastings ratio and accept or reject
 * 4) If step is accepted then reweight/update
 * 5) Sample from the Wishart for the variance matrix.
 * 6) If necessary, sample from the scale factor distribution (not done yet).
*/

/**
 * Create a Markov Chain Monte Carlo sample from a fitted generalized
 * linear mixed model
 *
 * @param GSpt External pointer to a GlmerStruct
 * @param savebp Logical indicator of whether or not to save the
 *   random effects in the MCMC sample
 * @param nsampp number of samples to generate
 * @param transp pointer to an logical scalar indicating if the
 * variance components should be transformed.
 *
 * @return a matrix
 */
SEXP
glmer_MCMCsamp(SEXP GSpt, SEXP savebp, SEXP nsampp, SEXP transp, SEXP verbosep)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSpt);
    SEXP ans, x = GS->mer;
    SEXP Omega = GET_SLOT(x, lme4_OmegaSym),
	Omegacp = PROTECT(duplicate(Omega));
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	i, j, nf = LENGTH(Omega), nsamp = asInteger(nsampp),
	p = LENGTH(GET_SLOT(x, lme4_rXySym)),
	q = LENGTH(GET_SLOT(x, lme4_rZySym)),
	saveb = asLogical(savebp),
	trans = asLogical(transp),
	verbose = asLogical(verbosep);
    double
	*RXX = REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym), lme4_xSym)),
	*RZX = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXSym), lme4_xSym)),
	*bhat = REAL(GET_SLOT(x, lme4_ranefSym)),
	*betahat = REAL(GET_SLOT(x, lme4_fixefSym)),
	*bcur = Calloc(q, double), *betacur = Calloc(p, double), /* current */
	*bnew = Calloc(q, double), *betanew = Calloc(p, double), /* proposed */
	*ansp, MHratio;
    int nrbase = p + 1 + coef_length(nf, nc); /* rows always included */
    int nrtot = nrbase + (saveb ? q : 0);

    if (nsamp <= 0) nsamp = 1;
    ans = PROTECT(allocMatrix(REALSXP, nrtot, nsamp));
    ansp = REAL(ans);
    for (i = 0; i < nrtot * nsamp; i++) ansp[i] = NA_REAL;
    GetRNGstate();
    /* initialize current values of b and beta to estimates */
    Memcpy(betacur, REAL(GET_SLOT(x, lme4_fixefSym)), p);
    Memcpy(bcur, REAL(GET_SLOT(x, lme4_ranefSym)), q);
    if(verbose) Rprintf("%12s\n", "MHratio");
    for (i = 0; i < nsamp; i++) {
	double odev, ndev, ogaus, ngaus, *col = ansp + i * nrtot;
				/* update bhat */
/* 	internal_bhat(GS, betahat, (double *) NULL); */
				/* generate proposal for b and beta */
				/* save Gaussian deviance at proposal */
	ngaus = internal_betab_update(p, q, 1, L, RZX, RXX,
				      betahat, bhat, betanew, bnew);
				/* Gaussian deviance at current */
	ogaus = internal_Gaussian_deviance(p, q, L, RZX, RXX, betahat, bhat,
					 betacur, bcur);
				/*  glmm deviance at proposal */
	ndev = glmm_deviance(GS, betanew, bnew);
				/*  glmm deviance at current */
	odev = glmm_deviance(GS, betacur, bcur);
				/* evaluate Metropolis-Hastings ratio */
	MHratio = exp((ngaus - ogaus - ndev + odev)/2);
	if(verbose) Rprintf("%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f\n",
			    ogaus, ngaus, odev, ndev, ngaus - ogaus, odev - ndev,
			    MHratio);
	/* Accept proposal with probability min(MHratio, 1) */
	if (unif_rand() < MHratio) {
	    Memcpy(betacur, betanew, p);
	    Memcpy(bcur, bnew, q);
	}
				/* Store beta */
	for (j = 0; j < p; j++) col[j] = betacur[j];
				/* Optionally store b */
	if (saveb) for (j = 0; j < q; j++) col[nrbase + j] = bcur[j];
				/* Update and store Omega */
	internal_Omega_update(Omega, bcur, 1, nf, nc, Gp, col + p, trans);
	internal_mer_refactor(x);
	mer_secondary(x);

	col[nrbase - 1] = glmm_deviance(GS, betacur, bcur); /* store deviance */
    }
    PutRNGstate();
    Free(bcur); Free(betacur); Free(betanew); Free(bnew);
				/* Restore original Omega */
    SET_SLOT(x, lme4_OmegaSym, Omegacp);
    internal_bhat(GS, betahat, (double *) NULL);
    internal_mer_refactor(x);

    Free(L);
    UNPROTECT(2);
    return ans;
}

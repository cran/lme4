#include "lme4_utils.h"

SEXP attr_hidden
alloc_dgeMatrix(int m, int n, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2));

    dims[0] = m; dims[1] = n;
    ALLOC_SLOT(ans, lme4_xSym, REALSXP, m * n);
    dn = ALLOC_SLOT(ans, lme4_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

SEXP attr_hidden
alloc_dpoMatrix(int n, char *uplo, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dpoMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, lme4_xSym, REALSXP, n * n);
    SET_SLOT(ans, lme4_uploSym, mkString(uplo));
    dn = ALLOC_SLOT(ans, lme4_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

SEXP attr_hidden
alloc_dtrMatrix(int n, char *uplo, char *diag, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtrMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, lme4_xSym, REALSXP, n * n);
    SET_SLOT(ans, lme4_uploSym, mkString(uplo));
    SET_SLOT(ans, lme4_diagSym, mkString(diag));
    dn = ALLOC_SLOT(ans, lme4_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

SEXP attr_hidden
alloc_dsCMatrix(int n, int nz, char *uplo, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsCMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, lme4_xSym, REALSXP, nz);
    ALLOC_SLOT(ans, lme4_iSym, INTSXP, nz);
    ALLOC_SLOT(ans, lme4_pSym, INTSXP, n + 1);
    SET_SLOT(ans, lme4_uploSym, mkString(uplo));
    dn = ALLOC_SLOT(ans, lme4_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

SEXP attr_hidden
alloc_dgCMatrix(int m, int n, int nz, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2));

    dims[0] = m; dims[1] = n;
    ALLOC_SLOT(ans, lme4_xSym, REALSXP, nz);
    ALLOC_SLOT(ans, lme4_iSym, INTSXP, nz);
    ALLOC_SLOT(ans, lme4_pSym, INTSXP, n + 1);
    dn = ALLOC_SLOT(ans, lme4_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
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
SEXP attr_hidden
alloc3Darray(SEXPTYPE mode, int nrow, int ncol, int nface)
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
 * Re-evaluate the factorization of the components of Omega and the
 * factorization of the crossproduct matrix - usually after an update.
 *
 * @param x pointer to an mer object
 */
/* NOTE: The factorization of the Omega components should be made
 * part of mer_factor.  Add an argument to mer_factor called force. */
void attr_hidden
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


/**
 * Calculate fitted values for the current fixed and random effects.
 *
 * @param x pointer to an mer object
 * @param initial initial values (i.e. an offset) or (double *) NULL
 * @param val array to hold the fitted values
 *
 * @return pointer to a numeric array of fitted values
 */
double attr_hidden
*internal_mer_fitted(SEXP x, const double initial[], double val[])
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
 * Set the coefficient vector and perform a factorization
 *
 * @param x pointer to an mer object
 * @param cc vector of coefficients to assign
 * @param ptyp indicator of the parametrization being used
 */
void attr_hidden internal_mer_coefGets(SEXP x, const double cc[], int ptyp)
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

/**
 * Inflate Z'Z to Z'Z+Omega and factor. Form RZX and rZy and update
 * the status flags.
 *
 * @param x pointer to an mer object.
 */
void attr_hidden internal_mer_Zfactor(SEXP x, cholmod_factor *L)
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
int attr_hidden internal_mer_Xfactor(SEXP x)
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
double attr_hidden
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

/**
 * Update the ranef slot assuming that the fixef, rZy, RZX and L slots
 * are up to date.
 *
 * @param x Pointer to an mer object
 *
 * @return pointer to the updated modes of the random effects
 */
double attr_hidden *internal_mer_ranef(SEXP x)
{
    SEXP ranef = GET_SLOT(x, lme4_ranefSym);
    int *status = INTEGER(GET_SLOT(x, lme4_statusSym));
    double *ans = REAL(ranef);

    if (!status[0]) {
	error("Applying internal_mer_ranef without factoring");
	return (double*)NULL;	/* -Wall */
    }
    if (status[0] < 2) {
	SEXP fixef = GET_SLOT(x, lme4_fixefSym);
	int ione = 1, p = LENGTH(fixef), q = LENGTH(ranef);
	cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
	cholmod_dense *td1, *td2, *chranef = M_as_cholmod_dense(ranef);
	double *RZX = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXSym), lme4_xSym)),
	    m1[] = {-1,0}, one[] = {1,0};

	Memcpy(ans, REAL(GET_SLOT(x, lme4_rZySym)), q);
	F77_CALL(dgemv)("N", &q, &p, m1, RZX, &q, REAL(fixef), &ione,
			one, ans, &ione);
	td1 = M_cholmod_solve(CHOLMOD_Lt, L, chranef, &c);
	td2 = M_cholmod_solve(CHOLMOD_Pt, L, td1, &c);
	Free(chranef); M_cholmod_free_dense(&td1, &c);
	Memcpy(ans, (double *)(td2->x), q);
	M_cholmod_free_dense(&td2, &c);
	status[0] = 2;
	Free(L);
    }
    return ans;
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
void attr_hidden
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
	std_rWishart_factor((double) (nlev - nci + 1), nci, 1, wfac);

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
void attr_hidden internal_ECMEsteps(SEXP x, int nEM, int verb)
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



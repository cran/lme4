
#if 0

/* Gauss-Hermite Quadrature x positions and weights */
static const double
    GHQ_x1[1] = {0},
    GHQ_w1[1] = {1},
    GHQ_x2[1] = {1},
    GHQ_w2[1] = {0.5},
    GHQ_x3[2] = {1.7320507779261, 0},
    GHQ_w3[2] = {0.166666666666667, 0.666666666666667},
    GHQ_x4[2] = {2.3344141783872, 0.74196377160456},
    GHQ_w4[2] = {0.0458758533899086, 0.454124131589555},
    GHQ_x5[3] = {2.85696996497785, 1.35562615677371, 0},
    GHQ_w5[3] = {0.0112574109895360, 0.222075915334214,
		 0.533333317311434},
    GHQ_x6[3] = {3.32425737665988, 1.88917584542184,
		 0.61670657963811},
    GHQ_w6[3] = {0.00255578432527774, 0.0886157433798025,
		 0.408828457274383},
    GHQ_x7[4] = {3.7504396535397, 2.36675937022918,
		 1.15440537498316, 0},
    GHQ_w7[4] = {0.000548268839501628, 0.0307571230436095,
		 0.240123171391455, 0.457142843409801},
    GHQ_x8[4] = {4.14454711519499, 2.80248581332504,
		 1.63651901442728, 0.539079802125417},
    GHQ_w8[4] = {0.000112614534992306, 0.0096352198313359,
		 0.117239904139746, 0.373012246473389},
    GHQ_x9[5] = {4.51274578616743, 3.20542894799789,
		 2.07684794313409, 1.02325564627686, 0},
    GHQ_w9[5] = {2.23458433364535e-05, 0.00278914123744297,
		 0.0499164052656755, 0.244097495561989,
		 0.406349194142045},
    GHQ_x10[5] = {4.85946274516615, 3.58182342225163,
		  2.48432579912153, 1.46598906930182,
		  0.484935699216176},
    GHQ_w10[5] = {4.31065250122166e-06, 0.000758070911538954,
		  0.0191115799266379, 0.135483698910192,
		  0.344642324578594},
    GHQ_x11[6] = {5.18800113558601, 3.93616653976536,
		  2.86512311160915, 1.87603498804787,
		  0.928868981484148, 0},
    GHQ_w11[6] = {8.12184954622583e-07, 0.000195671924393029,
		  0.0067202850336527, 0.066138744084179,
		  0.242240292596812, 0.36940835831095};

static const double
    *GHQ_x[12] = {(double *) NULL, GHQ_x1, GHQ_x2, GHQ_x3, GHQ_x4,
		  GHQ_x5, GHQ_x6, GHQ_x7, GHQ_x8, GHQ_x9, GHQ_x10,
		  GHQ_x11},
    *GHQ_w[12] = {(double *) NULL, GHQ_w1, GHQ_w2, GHQ_w3, GHQ_w4,
		  GHQ_w5, GHQ_w6, GHQ_w7, GHQ_w8, GHQ_w9, GHQ_w10,
		  GHQ_w11};

/**
 * Evaluate the conditional deviance for a given set of fixed effects.
 *
 * @param GS Pointer to a GlmerStruct
 * @param fixed value of the fixed effects
 *
 * @return conditional deviance
 */
static double
fixed_effects_deviance(GlmerStruct GS, const double fixed[])
{
    SEXP devs;
    int i, ione = 1;
    double ans, one = 1, zero = 0;

    F77_CALL(dgemv)("N", &(GS->n), &(GS->p), &one,
		    GS->X, &(GS->n), fixed,
		    &ione, &zero, REAL(GS->eta), &ione);
				/* add in random effects and offset */
    vecIncrement(REAL(GS->eta), GS->off, GS->n);
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    devs = PROTECT(eval_check(GS->dev_resids, GS->rho, REALSXP, GS->n));
    for (i = 0, ans = 0; i < GS->n; i++) ans += REAL(devs)[i];
    UNPROTECT(1);
    return ans;
}


/**
 * Evaluate the quadratic form (x-mn)'A'A(x-mn) from the multivariate
 * normal kernel.
 *
 * @param n dimension of random variate
 * @param mn mean vector
 * @param a upper Cholesky factor of the precision matrix
 * @param lda leading dimension of A
 * @param x vector at which to evaluate the kernel
 *
 * @return value of the normal kernel
 */
static double
normal_kernel(int n, const double mn[],
	      const double a[], int lda, const double x[])
{
    int i, ione = 1;
    double *tmp = Calloc(n, double), ans;

    for (i = 0; i < n; i++) tmp[i] = x[i] - mn[i];
    F77_CALL(dtrmv)("U", "N", "N", &n, a, &lda, tmp, &ione);
    for (i = 0, ans = 0; i < n; i++) ans += tmp[i] * tmp[i];
    Free(tmp);
    return ans;
}

/* FIXME: Separate the calculation of the offset random effects from
 * the calculation of the deviance contributions. */

/**
 * Determine the deviance components associated with each of the
 * levels of a grouping factor at the conditional modes or a value
 * offset from the conditional modes by delb.
 *
 * @param GS pointer to a GlmerStruct
 * @param b conditional modes of the random effects
 * @param Gp group pointers
 * @param nc number of columns in the model matrix for the kth
 * grouping factor
 * @param k index (0-based) of the grouping factor
 * @param delb vector of length nc giving the changes in the
 * orthonormalized random effects
 * @param OmgFac Cholesky factor of the inverse of the penalty matrix
 * for this grouping factor
 * @param bVfac 3-dimensional array holding the factors of the
 * conditional variance-covariance matrix of the random effects
FIXME: This is wrong.  It is bVar[[i]] not bVfac that is being passed.
This only affects the AGQ method.
 * @param devcmp array to hold the deviance components
 *
 * @return devcmp
 */
static double*
rel_dev_1(GlmerStruct GS, const double b[], int nlev, int nc, int k,
	  const double delb[], const double OmgFac[],
	  const double bVfac[], double devcmp[])
{
    SEXP devs;
    int *fv = INTEGER(VECTOR_ELT(GET_SLOT(GS->mer, lme4_flistSym), k)),
	i, j;
    double *bcp = (double *) NULL;

    AZERO(devcmp, nlev);
    if (delb) {
	int ione = 1, ntot = nlev * nc;
	double sumsq = 0;
				/* copy the contents of b */
	bcp = Memcpy(Calloc(ntot, double), b, ntot);
	if (nc == 1) {
	    sumsq = delb[0] * delb[0];
	    for (i = 0; i < nlev; i++) b[i] += delb[0] * bVfac[i];
	} else {
	    int ncsq = nc * nc;
	    double *tmp = Calloc(nc, double);
	    for (i = 0; i < nlev; i++) {
		Memcpy(tmp, delb, nc);
		F77_CALL(dtrmv)("U", "N", "N", &nc, &(bVfac[i * ncsq]),
				&nc, tmp, &ione);
		for (j = 0; j < nc; j++) b[i + j * nc] = tmp[j];
	    }
				/* sum of squares of delb */
	    for (j = 0; j < nc; j++) sumsq += delb[j] * delb[j];
	}
	for (i = 0; i < nlev; i++) devcmp[i] = -sumsq;
    }
    internal_mer_fitted(GS->mer, GS->offset, REAL(GS->eta));
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    devs = PROTECT(eval_check(GS->dev_resids, GS->rho, REALSXP, GS->n));
    for (i = 0; i < GS->n; i++)
	devcmp[fv[i] - 1] += REAL(devs)[i];
    UNPROTECT(1);
    if (nc == 1) {
	for (i = 0; i < nlev; i++) {
	    double tmp = *OmgFac * b[i];
	    devcmp[i] += tmp * tmp;
	}
    } else {
	double *tmp = Calloc(nc, double);
	int ione = 1;

	for (i = 0; i < nlev; i++) {
	    for (j = 0; j < nc; j++) tmp[j] = b[i + j * nlev];
	    F77_CALL(dtrmv)("U", "N", "N", &nc, OmgFac, &nc,
			    tmp, &ione);
	    for (j = 0; j < nc; j++)
		devcmp[i] += tmp[j] * tmp[j];
	}
    }
    if (delb) {
	Memcpy(b, bcp, ntot);
	Free(bcp);
    }
    return devcmp;
}


/**
 * Determine the conditional modes and the conditional variance of the
 * fixed effects given the data and the current random effects.
 * Create a Metropolis-Hasting proposal step from the multivariate
 * normal density, determine the acceptance probability and return the
 * current value or the proposed value.
 *
 * @param GSp pointer to a GlmerStruct
 * @param b list of random effects
 * @param fixed current value of the fixed effects
 *
 * @return updated value of the fixed effects
 */
SEXP glmer_fixed_update(SEXP GSp, SEXP b, SEXP fixed)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);

    if (!isReal(fixed) || LENGTH(fixed) != GS->p)
	error(_("%s must be a %s of length %d"), "fixed",
		"numeric vector", GS->p);
    GetRNGstate();
    internal_glmer_fixef_update(GS, b, REAL(fixed));
    PutRNGstate();
    return fixed;
}

SEXP glmer_bhat(SEXP pars, SEXP GSp)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);

    if (!isReal(pars) || LENGTH(pars) != GS->npar)
	error(_("`%s' must be a numeric vector of length %d"),
	      "pars", GS->npar);
    if (!internal_bhat(GS, REAL(pars), REAL(pars) + (GS->p)))
	warning(_("internal_bhat did not converge"));
    return R_NilValue;
}



/**
 * Determine the conditional modes and the conditional variance of the
 * fixed effects given the data and the current random effects.
 * Create a Metropolis-Hasting proposal step from the multivariate
 * normal density, determine the acceptance probability and, if the
 * step is to be accepted, overwrite the contents of fixed with the
 * new contents.
 *
 * @param GS a GlmerStruct
 * @param b list of random effects
 * @param fixed current value of the fixed effects
 *
 * @return updated value of the fixed effects
 */
static double *
internal_glmer_fixef_update(GlmerStruct GS, SEXP b,
			    double fixed[])
{
    SEXP dmu_deta, var;
    int i, ione = 1, it, j, lwork = -1;
    double *ans = Calloc(GS->p, double), /* proposal point */
	*md = Calloc(GS->p, double), /* conditional modes */
	*w = Calloc(GS->n, double), *work,
	*wtd = Calloc(GS->n * GS->p, double),
	*z = Calloc(GS->n, double),
	crit, devr, one = 1, tmp, zero = 0;

    if (!isNewList(b) || LENGTH(b) != GS->nf)
	error(_("%s must be a %s of length %d"), "b", "list", GS->nf);
    for (i = 0; i < GS->nf; i++) {
	SEXP bi = VECTOR_ELT(b, i);
	if (!isReal(bi) || !isMatrix(bi))
	    error(_("b[[%d]] must be a numeric matrix"), i);
    }
    AZERO(z, GS->n);		/* -Wall */
    Memcpy(md, fixed, GS->p);
				/* calculate optimal size of work array */
    F77_CALL(dgels)("N", &(GS->n), &(GS->p), &ione, wtd, &(GS->n),
		    z,  &(GS->n), &tmp, &lwork, &j);
    if (j)			/* shouldn't happen */
	error(_("%s returned error code %d"), "dgels", j);
    lwork = (int) tmp;
    work = Calloc(lwork, double);

    AZERO(GS->off, GS->n); /* fitted values from random effects */
/*     fitted_ranef(GET_SLOT(GS->mer, lme4_flistSym), GS->unwtd, b, */
/* 		 INTEGER(GET_SLOT(GS->mer, lme4_ncSym)), GS->off); */
    for (i = 0; i < GS->n; i++)
	(GS->etaold)[i] = ((GS->off)[i] += (GS->offset)[i]);

    for (it = 0, crit = GS->tol + 1;
	 it < GS->maxiter && crit > GS->tol; it++) {
				/* fitted values from current beta */
	F77_CALL(dgemv)("N", &(GS->n), &(GS->p), &one,
			GS->X, &(GS->n), md,
			&ione, &zero, REAL(GS->eta), &ione);
				/* add in random effects and offset */
	vecIncrement(REAL(GS->eta), (GS->off), GS->n);
				/* check for convergence */
	crit = conv_crit(GS->etaold, REAL(GS->eta), GS->n);
				/* obtain mu, dmu_deta, var */
	eval_check_store(GS->linkinv, GS->rho, GS->mu);
	dmu_deta = PROTECT(eval_check(GS->mu_eta, GS->rho,
				      REALSXP, GS->n));
	var = PROTECT(eval_check(GS->var, GS->rho, REALSXP, GS->n));
				/* calculate weights and working residual */
	for (i = 0; i < GS->n; i++) {
	    w[i] = GS->wts[i] * REAL(dmu_deta)[i]/sqrt(REAL(var)[i]);
	    z[i] = w[i] * (REAL(GS->eta)[i] - (GS->off)[i] +
			   ((GS->y)[i] - REAL(GS->mu)[i]) /
			   REAL(dmu_deta)[i]);
	}
	UNPROTECT(2);
				/* weighted copy of the model matrix */
	for (j = 0; j < GS->p; j++)
	    for (i = 0; i < GS->n; i++)
		wtd[i + j * GS->n] = GS->X[i + j * GS->n] * w[i];
				/* weighted least squares solution */
	F77_CALL(dgels)("N", &(GS->n), &(GS->p), &ione, wtd, &(GS->n),
			z, &(GS->n), work, &lwork, &j);
	if (j) error(_("%s returned error code %d"), "dgels", j);
	Memcpy(md, z, GS->p);
    }
				/* wtd contains the Cholesky factor of
				 * the precision matrix */
    devr = normal_kernel(GS->p, md, wtd, GS->n, fixed);
    devr -= fixed_effects_deviance(GS, fixed);
    for (i = 0; i < GS->p; i++) {
	double var = norm_rand();
	ans[i] = var;
	devr -= var * var;
    }
    F77_CALL(dtrsv)("U", "N", "N", &(GS->p), wtd, &(GS->n), ans, &ione);
    for (i = 0; i < GS->p; i++) ans[i] += md[i];
    devr += fixed_effects_deviance(GS, ans);
    crit = exp(-0.5 * devr);	/* acceptance probability */
    tmp = unif_rand();
    if (asLogical(internal_getElement(GS->cv, "msVerbose"))) {
	Rprintf("%5.3f: ", crit);
	for (j = 0; j < GS->p; j++) Rprintf("%#10g ", ans[j]);
	Rprintf("\n");
    }
    if (tmp < crit) Memcpy(fixed, ans, GS->p);
    Free(ans); Free(md); Free(w);
    Free(work); Free(wtd); Free(z);
    return fixed;
}

static void
internal_glmer_ranef_update(GlmerStruct GS, SEXP b)
{
    SEXP bhat, bprop = PROTECT(duplicate(b)),
	bVar = GET_SLOT(GS->mer, lme4_bVarSym),
	Omega = GET_SLOT(GS->mer, lme4_OmegaSym);
    int i, ione = 1, j, k;
    double *b = REAL(GET_SLOT(GS->mer, lme4_ranefSym);
    double devr, one = 1;

    bhat = R_NilValue;
				/* subtract deviance at b */
    devr = -random_effects_deviance(GS, b);
     for (i = 0; i < GS->nf; i++) {
 	SEXP Bi = VECTOR_ELT(b, i);
 	int *dims = INTEGER(getAttrib(Bi, R_DimSymbol));
 	int nlev = dims[0], nci = dims[1];
 	int ncisqr = nci * nci, ntot = nlev * nci;
 	double *bcopy = Calloc(ntot, double),
 	    *bi = REAL(Bi),
 	    *bhati = REAL(VECTOR_ELT(bhat, i)),
 	    *bpropi = REAL(VECTOR_ELT(bprop, i)),
 	    *bVari = REAL(VECTOR_ELT(bVar, i)),
 	    *chol = Calloc(ncisqr, double),
 	    *delta = Calloc(nci, double),
 	    *omgfac = Memcpy(Calloc(ncisqr, double),
 			     REAL(VECTOR_ELT(Omega, i)),
 			     ncisqr);

 				/* subtract quadratic form in */
 				/* Omega[[i]] at b  */
 	F77_CALL(dpotrf)("U", &nci, omgfac, &nci, &j);
 	if (j)
 	    error(_("Leading %d minor of Omega[[%d]] not positive definite"),
                       j, i + 1);
 	Memcpy(bcopy, bi, ntot);
 	F77_CALL(dtrmm)("R", "U", "T", "N", &nlev, &nci, &one,
 			omgfac, &nci, bcopy, &nlev);
 	for (k = 0; k < ntot; k++) devr -= bcopy[k] * bcopy[k];
 				/* form bprop and proposal density */
 	for (k = 0; k < nlev; k++) {
 				/* proposal density at b */
 	    for (j = 0; j < nci; j++)
 		delta[j] = bi[k + j * nlev] - bhati[k + j * nlev];
 	    Memcpy(chol, &(bVari[k * ncisqr]), ncisqr);
 	    F77_CALL(dpotrf)("U", &nci, chol, &nci, &j);
 	    if (j)
 		error(_("Leading %d minor of bVar[[%d]][,,%d] not positive definite"),
                       j, i + 1, k + 1);
 	    F77_CALL(dtrsv)("U", "T", "N", &nci, chol, &nci,
 			    delta, &ione);
 	    for (j = 0; j < nci; j++) {
 		double nrm = norm_rand(); /* proposal deviate */
 		devr += delta[j] * delta[j] - nrm * nrm;
 		delta[j] = nrm;
 	    }
 				/* scale by Cholesky inverse */
 	    F77_CALL(dtrmv)("U", "T", "N", &nci, chol, &nci,
 			    delta, &ione);
 				/* add mean */
 	    for (j = 0; j < nci; j++)
 		bpropi[k + j * nlev] = bhati[k + j * nlev] + delta[j];
 	}
 				/* add quadratic form in */
 				/* Omega[[i]] at bprop  */
 	Memcpy(bcopy, bpropi, ntot);
 	F77_CALL(dtrmm)("R", "U", "T", "N", &nlev, &nci, &one,
 			omgfac, &nci, bcopy, &nlev);
 	for (k = 0; k < ntot; k++) devr += bcopy[k] * bcopy[k];

 	Free(bcopy); Free(chol); Free(delta); Free(omgfac);
     }
 				/* add deviance at bprop */
/*     devr += random_effects_deviance(GS, bprop); */

    if (unif_rand() < exp(-0.5 * devr))
	for (i = 0; i < GS->nf; i++) { /* copy each face of b */
	    SEXP Bi = VECTOR_ELT(b, i);
	    int *dims = INTEGER(getAttrib(Bi, R_DimSymbol));

	    Memcpy(REAL(Bi), REAL(VECTOR_ELT(bprop, i)),
		   dims[0] * dims[1]);
	}

    if (asLogical(internal_getElement(GS->cv, "msVerbose"))) {
	double *b0 = REAL(VECTOR_ELT(bprop, 0));
	Rprintf("%5.3f:", exp(-0.5 * devr));
	for (k = 0; k < 5; k++) Rprintf("%#10g ", b0[k]);
	Rprintf("\n");
    }

    UNPROTECT(2);
}


/**
 * Compute the approximation to the deviance using adaptive
 * Gauss-Hermite quadrature (AGQ).  When nAGQ == 1 this is the Laplace
 * approximation.
 *
 * @param pars pointer to a numeric vector of parameters
 * @param GSp pointer to a GlmerStruct object
 * @param nAGQp pointer to a scalar integer representing the number of
 * points in AGQ to use
 *
 * @return the approximation to the deviance as computed using AGQ
 */
SEXP glmer_devAGQ(SEXP pars, SEXP GSp, SEXP nAGQp)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    SEXP Omega = GET_SLOT(GS->mer, lme4_OmegaSym),
	bVar = GET_SLOT(GS->mer, lme4_bVarSym);
    int i, j, k, nAGQ = asInteger(nAGQp);
    int n2 = (nAGQ + 1)/2,
	*Gp = INTEGER(GET_SLOT(GS->mer, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(GS->mer, lme4_ncSym));
    double *f0, LaplaceDev = 0, AGQadjst = 0,
	*bhat = REAL(GET_SLOT(GS->mer, lme4_ranefSym));

    if (!isReal(pars) || LENGTH(pars) != GS->npar)
	error(_("`%s' must be a numeric vector of length %d"),
	      "pars", GS->npar);
    if (GS->nf > 1 && nAGQ > 1) {
	warning(_("AGQ not available for multiple grouping factors - using Laplace"));
	nAGQ = 1;
    }
    if (!internal_bhat(GS, REAL(pars), REAL(pars) + (GS->p)))
	return ScalarReal(DBL_MAX);

    for (i = 0; i < GS->nf; i++) {
	int nci = nc[i];
	int ncip1 = nci + 1, ncisqr = nci * nci,
	    nlev = (Gp[i + 1] - Gp[i])/nci;
	double *omgf = REAL(GET_SLOT(M_dpoMatrix_chol(VECTOR_ELT(Omega, i)), lme4_xSym)),
	    *bVi = Memcpy(Calloc(ncisqr * nlev, double),
			   REAL(VECTOR_ELT(bVar, i)), ncisqr * nlev);

        for (j = 0; j < nci; j++) { /* nlev * logDet(Omega_i) */
            LaplaceDev += 2 * nlev * log(omgf[j * ncip1]);
        }
        for (k = 0; k < nlev; k++) {
	    double *bVik = bVi + k * ncisqr;
            F77_CALL(dpotrf)("U", &nci, bVik, &nci, &j);
            if (j)
                error(_("Leading %d minor of bVar[[%d]][,,%d] not positive definite"),
                      j, i + 1, k + 1);
            for (j = 0; j < nci; j++) LaplaceDev -= 2 * log(bVik[j * ncip1]);
        }

	f0 = Calloc(nlev, double);
	rel_dev_1(GS, bhat, nlev, nci, i, (double *) NULL,
		  omgf, bVi, f0);
	for (k = 0; k < nlev; k++) LaplaceDev += f0[k];
	if (nAGQ > 1) {
	    double *fx = Calloc(nlev, double),
		*rellik = Calloc(nlev, double),
		*delb = Calloc(nci, double);

	    if (nci > 1) error(_("code not yet written"));
	    AZERO(rellik, nlev);	/* zero accumulator */
	    for (k = 0; k < n2; k++) {
		delb[0] = GHQ_x[nAGQ][k];
		if (delb[0]) {
		    rel_dev_1(GS, bhat, nlev, nci, i, delb,
			      omgf, bVi, fx);
		    for (j = 0; j < nlev; j++) {
			rellik[j] += GHQ_w[nAGQ][k] *
			    exp(-(fx[j] - f0[j])/2);
		    }
		    delb[0] *= -1;
		    rel_dev_1(GS, bhat, nlev, nci, i, delb,
			      omgf, bVi, fx);
		    for (j = 0; j < nlev; j++) {
			rellik[j] += GHQ_w[nAGQ][k] *
			    exp(-(fx[j] - f0[j])/2);
		    }
		} else {
		    for (j = 0; j < nlev; j++)
			rellik[j] += GHQ_w[nAGQ][k];
		}
	    }
	    for (j = 0; j < nlev; j++)
		AGQadjst -= 2 * log(rellik[j]);
	    Free(fx); Free(rellik);
	}
	Free(f0); Free(bVi);
    }
    return ScalarReal(LaplaceDev + AGQadjst);
}


/**
 * Determine the conditional modes and the conditional variance of the
 * random effects given the data and the current fixed effects and
 * variance components.
 *
 * Create a Metropolis-Hasting proposal step from a multivariate
 * normal density centered at bhat with variance-covariance matrix
 * from bVar, determine the acceptance probability and return the
 * current value or the proposed value.
 *
 * @param GSp pointer to a GlmerStruct
 * @param fixed pointer to a numeric vector of the fixed effects
 * @param varc pointer to a numeric vector of the variance components
 * @param varc pointer to current values of b
 *
 * @return updated b from the Metropolis-Hastings step
 */
SEXP glmer_ranef_update(SEXP GSp, SEXP fixed, SEXP varc, SEXP b)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    int nvarc = GS->npar - GS->p;

    if (!isReal(fixed) || LENGTH(fixed) != GS->p)
	error(_("`%s' must be a numeric vector of length %d"),
	      "fixed", GS->p);
    if (INTEGER(GET_SLOT(GS->mer, lme4_ncSym))[GS->nf] > 0)
	error(_("the mer object must be set to skip fixed effects"));
    if (!isReal(varc) || LENGTH(varc) != nvarc)
	error(_("`%s' must be a numeric vector of length %d"),
	      "varc", nvarc);

    GetRNGstate();
    /* Don't check for convergence failure after internal_bhat.
     * It is determining the mean of the proposal density and
     * does not need exact convergence. */
    internal_bhat(GS, REAL(fixed), REAL(varc));
    internal_glmer_ranef_update(GS, b);
    PutRNGstate();

    UNPROTECT(1);
    return b;
}


#endif

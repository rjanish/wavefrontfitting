/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nr.h.  Do not confuse this file with the same-named
   file nr.h that is supplied in the 'misc' subdirectory.
   *That* file is the one from the book, and contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#ifndef _NR_NEWT_H_
#define _NR_NEWT_H_

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {double r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

#ifndef _ARITHCODE_DECLARE_T_
typedef struct {
	unsigned long *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;
#define _ARITHCODE_DECLARE_T_
#endif /* _ARITHCODE_DECLARE_T_ */

#ifndef _HUFFCODE_DECLARE_T_
typedef struct {
	unsigned long *icod,*ncod,*left,*right,nch,nodemax;
} huffcode;
#define _HUFFCODE_DECLARE_T_
#endif /* _HUFFCODE_DECLARE_T_ */

#include <stdio.h>

// nonlinear system solver routines
void newt(double x[], int n, int *check, 
    double parameters[], double x_image, double y_image,
    void (*vecfunc)(int, double [], double [], double [], double, double));
void fdjac(int n, double x[], double fvec[], double **df,
    double parameters, double x_im, double y_im,
	void (*vecfunc)(int, double [], double [], double [], double, double));
void lubksb(double **a, int n, int *indx, double b[]);
void ludcmp(double **a, int n, int *indx, double *d);
void lnsrch(int n, double xold[], double fold, double g[], double p[], 
            double x[], double *f, double stpmax, int *check,
            double image_params[], double xim, double yim, 
            double (*func)(double [], double [], double, double));
double nr_fmin(double x[], double params[], double x_image, double y_imag);
// Levenberg–Marquardt fitting routines
void mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **covar, double **alpha, double *chisq,
	void (*funcs)(double, double [], double *, double [], int), double *alamda, int verbose);
void mrqcof(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **alpha, double beta[], double *chisq,
	void (*funcs)(double, double [], double *, double [], int));
void gaussj(double **a, int n, double **b, int m);
void covsrt(double **covar, int ma, int ia[], int mfit);

#endif /* _NR_NEWT_H_ */

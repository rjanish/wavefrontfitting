#include <math.h>
#define NRANSI
#include "nr_newt.h"
#include "nrutil.h"
#define MAXITS 10
#define TOLF 1.0e-4
#define TOLMIN 1.0e-6
#define TOLX 1.0e-7
#define STPMX 100.0

int nn;
double *fvec;
void (*nrfuncv)(int n, double v[], double f[]);
#define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);\
	free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);return;}

void newt(double x[], int n, int *check, int *where, 
	void (*vecfunc)(int, double [], double []))
{
	void fdjac(int n, double x[], double fvec[], double **df,
		void (*vecfunc)(int, double [], double []));
	double nr_fmin(double x[]);
	void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
		 double *f, double stpmax, int *check, double (*func)(double []));
	void lubksb(double **a, int n, int *indx, double b[]);
	void ludcmp(double **a, int n, int *indx, double *d);
	int i,its,j,*indx;
	double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;

	indx=ivector(1,n);
	fjac=matrix(1,n,1,n);
	g=vector(1,n);
	p=vector(1,n);
	xold=vector(1,n);
	fvec=vector(1,n);
	nn=n;
	nrfuncv=vecfunc;
	f=nr_fmin(x);
	test=0.0;
	for (i=1;i<=n;i++)
		if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
	if (test < 0.01*TOLF) {
		*check=0;
		*where=1;
		FREERETURN
	}
	for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=MAXITS;its++) {
		fdjac(n,x,fvec,fjac,vecfunc);
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
			g[i]=sum;
		}
		for (i=1;i<=n;i++) xold[i]=x[i];
		fold=f;
		for (i=1;i<=n;i++) p[i] = -fvec[i];
		ludcmp(fjac,n,indx,&d);
		lubksb(fjac,n,indx,p);
		lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,nr_fmin);
		test=0.0;
		for (i=1;i<=n;i++)
			if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
		if (test < TOLF) {
			*check=0;
			*where=2;
			FREERETURN
		}
/*		if (*check) {
			test=0.0;
			den=FMAX(f,0.5*n);
			for (i=1;i<=n;i++) {
				temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
				if (temp > test) test=temp;
			}
			*check=(test < TOLMIN ? 1 : 0);
			*where = 3;
			FREERETURN
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) {
		    *where = 4;
		    FREERETURN
		} */
	}
	*where = 5;
	*check = 1;
	FREERETURN
}
#undef MAXITS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI


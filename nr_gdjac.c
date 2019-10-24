#include <math.h>
#define NRANSI
#include "nr_routines.h"
#include "nrutil.h"
#define EPS 1.0e-4

void fdjac(int n, double x[], double fvec[], double **df,
	double parameters[], double x_im, double y_im,
	void (*vecfunc)(int, double [], double [], double [], double, double))
{
	int i,j;
	double h,temp,*f;

	f=vector(1,n);
	for (j=1;j<=n;j++) {
		temp=x[j];
		h=EPS*fabs(temp);
		if (h == 0.0) h=EPS;
		x[j]=temp+h;
		h=x[j]-temp;
		(*vecfunc)(n,x,f, parameters, x_im, y_im);
		x[j]=temp;
		for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
	}
	free_vector(f,1,n);
}
#undef EPS
#undef NRANSI

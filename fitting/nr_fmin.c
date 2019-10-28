#define NRANSI
#include "nr_routines.h"
#include "nrutil.h"

extern int nn;
extern double *fvec;
extern void (*nrfuncv)(int n, double v[], double f[], double params[], 
                double x_image, double y_image);

double nr_fmin(double x[], double params[], double x_image, double y_image)
{
	int i;
	double sum;

	(*nrfuncv)(nn, x, fvec, params, x_image, y_image);
	for (sum=0.0,i=1;i<=nn;i++) sum += SQR(fvec[i]);
	return 0.5*sum;
}
#undef NRANSI

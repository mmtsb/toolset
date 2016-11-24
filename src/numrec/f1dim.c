#include <stdio.h>
#define NRANSI
#include "nrutil.h"

extern int ncom;
extern double *pcom,*xicom,(*nrfunc)(double []);

double f1dim(double x)
{
	int j;
	double f,*xt;

	xt=vector(1,ncom);
	for (j=1;j<=ncom;j++) {
          xt[j]=pcom[j]+x*xicom[j];
        }
	f=(*nrfunc)(xt);
	free_vector(xt,1,ncom);
	return f;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software *1(.|a. */

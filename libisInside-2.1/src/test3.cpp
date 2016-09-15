#include "isInside.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define npt   5
#define N  3000

#define eps ((double)(1e-30))
#define xmin  ((double)(0.0)-eps)
#define xmax  ((double)(0.0)+eps)
//#define ycst  ((double)(0.5))
#define ycst  ((double)(0.0))

int main(int narg, char **arg) {

  double **x;
  long i,j;
  signed_curve *c;
  FILE *f;
  double pt[2];

  x=(double **)malloc(npt*sizeof(double *));
  for (i=0;i<npt;i++) {
    x[i]=(double *)malloc(2*sizeof(double));
  }

  x[0][0]=(double)( 0.0);
  x[0][1]=(double)(-1.0);
  x[1][0]=(double)( 1.0);
  x[1][1]=(double)(-1.0);
  x[2][0]=(double)( 1.0);
  x[2][1]=(double)( 0.0);
  x[3][0]=(double)( 1.0);
  x[3][1]=(double)( 1.0);
  x[4][0]=(double)( 0.0);
  x[4][1]=(double)( 1.0);

  c=new signed_curve(npt,x);

  printf("AREA=%lf\n",c->area());

  f=fopen("test.dat","w");
  fprintf(f,"VARIABLES=\"x\"\"isIns\"\n");
  fprintf(f,"ZONE\n");
  pt[0]=xmin;
  pt[1]=ycst;
    for (i=0;i<=2*N+1;i++) {
      fprintf(f,"%30.25lf\t%30.25lf\n",
	      (1e10)*pt[0],
	      c->d_isInside(pt));
      pt[0]+=eps/((double)(N));
      pt[1]+=0.8/((double)(N));
    }

  fclose(f);

  return(0);
}

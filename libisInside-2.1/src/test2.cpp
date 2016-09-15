#include "isInside.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define npt 100
#define N   300

#define xmin -2.0
#define xmax  2.0
#define ymin -2.0
#define ymax  2.0

//#define xmin  0.0
//#define xmax  1.0
//#define ymin  0.0
//#define ymax  1.0

int main(int narg, char **arg) {

  double **x;
  long i,j;
  signed_curve *c;
  FILE *f;
  FILE *f2;
  double pt[2];

  x=(double **)malloc(npt*sizeof(double *));
  for (i=0;i<npt;i++) {
    x[i]=(double *)malloc(2*sizeof(double));
  }

  f2=fopen("boundary.dat","w");
  fprintf(f2,"%i\n",npt);
  for (i=0;i<npt;i++) {
    x[i][0]=cos(i*2*M_PI/npt);
    x[i][1]=sin(i*2*M_PI/npt);
    fprintf(f2,"%lf\t%lf\n",x[i][0],x[i][1]);
  }
  fclose(f2);

  c=new signed_curve("boundary.dat");

  printf("AREA=%lf\n",c->area());

  f=fopen("test.dat","w");
  fprintf(f,"VARIABLES=\"x\"\"y\"\"isIns\"\n");
  fprintf(f,"ZONE I=%i J=%i\n",N+1,N+1);
  for (j=0;j<=N;j++) {
    for (i=0;i<=N;i++) {
      pt[0]=xmin+(xmax-xmin)*((double)(i))/((double)(N));
      pt[1]=ymin+(ymax-ymin)*((double)(j))/((double)(N));
      fprintf(f,"%lf\t%lf\t%i\n",
	      pt[0],pt[1],
	      c->isInside(pt));
    }
  }

  fclose(f);

  return(0);
}

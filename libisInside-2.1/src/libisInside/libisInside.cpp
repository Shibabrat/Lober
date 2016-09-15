#include "isInside.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int signed_curve::rotation(void) {
  return(this->rot);
}

void signed_curve::make(long npt, double **x) {
  long i;
  this->npt=npt;
  pt[0]=(double *)malloc((npt+1)*sizeof(double));
  pt[1]=(double *)malloc((npt+1)*sizeof(double));
  for (i=0;i<npt;i++) {
    pt[0][i]=x[i][0];
    pt[1][i]=x[i][1];
  }
  pt[0][npt]=x[0][0];
  pt[1][npt]=x[0][1];
  myArea=area_compute();
  rot=1;
  if (myArea<(double)(0.0)) {
    myArea*=-1;
    rot=-1;
  }
}

void signed_curve::make(long npt, double *x) {
  long i;
  this->npt=npt;
  pt[0]=(double *)malloc((npt+1)*sizeof(double));
  pt[1]=(double *)malloc((npt+1)*sizeof(double));
  for (i=0;i<npt;i++) {
    pt[0][i]=x[2*i];
    pt[1][i]=x[2*i+1];
  }
  pt[0][npt]=x[0];
  pt[1][npt]=x[1];
  myArea=area_compute();
  rot=1;
  if (myArea<(double)(0.0)) {
    myArea*=-1;
    rot=-1;
  }
}

signed_curve::signed_curve(long npt, double *x) {
  this->make(npt,x);
}

signed_curve::signed_curve(long npt, double **x) {
  this->make(npt,x);
}

signed_curve::signed_curve(const char *filename) {
  FILE *f;
  double **q=NULL;
  long np=0;
  long i;
  double x,y;
  f=fopen(filename,"r");
  if (f==NULL) {
    printf("LIBISINSIDE::Cannot open file \"%s\"\n",filename);
  }
  // fscanf(f,"%i",&np);
  fscanf(f,"%ld",&np);
  //printf("n=%i\n",np);
  q=(double **)malloc(np*sizeof(double *));
  for (i=0;i<np;i++) {
    q[i]=(double *)malloc(2*sizeof(double));
    fscanf(f,"%lf",&x);
    fscanf(f,"%lf",&y);
    q[i][0]=x;
    q[i][1]=y;
    //printf("%45.40lf\t%45.40lf\n",x,y);
  }
  this->make(np,q);
  if (f!=NULL) {
    fclose(f);
  }
  for (i=0;i<np;i++) {
    free(q[i]);
  }
  if (q!=NULL) {
    free(q);
  }
}

signed_curve::~signed_curve() {
  free(pt[0]);
  free(pt[1]);
}

double signed_curve::area(void) {
  //printf("here\n");
  //printf("npt=%i\n",npt);
  return(myArea);
}


double signed_curve::area_compute(void) {
  double a=(double)(0.0);
  double x,y,dx,dy;
  double dA;
  long i;
  for (i=0;i<npt;i++) {
    x=.5*(pt[0][i]+pt[0][i+1]);
    y=.5*(pt[1][i]+pt[1][i+1]);
    dx=pt[0][i+1]-pt[0][i];
    dy=pt[1][i+1]-pt[1][i];
    dA=y*dx-x*dy;
    a+=0.5*dA;
  }
  return(a);
}

/***OLD VERSION***
double superlog(double x) {
	//filter for Math.h->log()
	if (x>(double)(0.0)) 
		return(log(x));
	else
		return((double)(0.0));
}

int signed_curve::isInside(double point[2]) {
  int isIns;
  double xint, x1, y1, x2, y2;
  double a,b,c,delta,factor;
  double normal;
  long i=0;

  double x,y;
  x=point[0];
  y=point[1];

  //printf("checking %f %f",x,y);
  isIns=1;
  xint=((double)(0.0));
  x2=pt[0][0];
  y2=pt[1][0];
  for (i=0;i<npt;i++) {
    x1=x2;
    y1=y2;
    x2=pt[0][i+1];
    y2=pt[1][i+1];

    a=(x1-x)*(x1-x)+(y1-y)*(y1-y);
    b=2.0*((x2-x1)*(x1-x)+(y2-y1)*(y1-y));
    c=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1);
    normal=fabs(a);
    normal=(fabs(b)>normal)?fabs(b):normal;
    normal=(fabs(c)>normal)?fabs(c):normal;
    // a.b.c are normalized 
    a/=normal;
    b/=normal;
    c/=normal;
    // delta : normalized 
    delta=4*a*c-b*b;
    factor=(x1-x)*(y2-y1)-(y1-y)*(x2-x1);
    if (delta<0.0) {
      delta=sqrt(-delta);
      xint+=(factor/normal)*(1.0/delta)*superlog(((b+2*c-delta)*(b+delta))/((b+2*c+delta)*(b-delta)));
      //sprintf(txt,"unexpected delta=%f\n",delta);if (rank==0) debugOut(txt);
      //exit(0);
    } else if (delta==0.0) {
      //xint-=2.0*(factor/normal)/(b+2*c);
      //xint+=2.0*(factor/normal)/b;
      //sprintf(txt,"unexpected delta=%f\n",delta);if (rank==0) debugOut(txt);
      //exit(0);
    } else {
      delta=sqrt(delta);
      xint+=(factor/normal)*(2.0/delta)*atan((b+2*c)/delta);
      xint-=(factor/normal)*(2.0/delta)*atan(b/delta);
    }
    
  }
  
  xint*=-rot;

  if (xint<3.1415926535)
    isIns=0;
  
  return(isIns);
}
***/

double signed_curve::d_isInside(double point[2]) {
  double xint, x1, y1, x2, y2;
  //double a,b,c;
  double delta,factor;
  double normal;
  long i=0;
  double uv, vv, uxv;

  double x,y;
  x=point[0];
  y=point[1];

  xint=((double)(0.0));
  x2=pt[0][0];
  y2=pt[1][0];
  for (i=0;i<npt;i++) {
    x1=x2;
    y1=y2;
    x2=pt[0][i+1];
    y2=pt[1][i+1];
    uv=(x2-x1)*(x1-x)+(y2-y1)*(y1-y);
    vv=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1);
    uxv=(x1-x)*(y2-y1)-(y1-y)*(x2-x1);
    if (vv!=((double)(0.0))) {
      normal=fabs(uv);
      normal=(fabs(vv)>normal)?fabs(vv):normal;
      normal=(fabs(uxv)>normal)?fabs(uxv):normal;
      //normalization
      uv/=normal;
      vv/=normal;
      uxv/=normal;
      if (uxv==(double)(0.0)) {
      //if (fabs(uxv)<(double)(1e-10)) {
	//u is parall to v. If the point is NOT
	//on a segment, not a big deal because the increment
	//to the integral is zero anyway.
	if ((x2-x)*(x1-x)+(y2-y)*(y1-y)<=(double)(0.0)) {
	  return(M_PI);
	}
      } else {
	xint+=atan((vv+uv)/uxv)-atan(uv/uxv);
      }    
    }
  }
  
  xint*=-rot;
  
  return(xint);
}

int signed_curve::isInside(double point[2]) {
  int isIns=1;
  double xint;
  
  xint=d_isInside(point);

  if (xint<M_PI/2.0) {
    isIns=-1;
  } else if (xint<3*M_PI/2.0) {
    isIns=0;
  }
  
  return(isIns);
}

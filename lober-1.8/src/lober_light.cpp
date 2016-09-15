#include "lober.h"
#include "config.h"
#include "inandout.h"
#include <isInside.h>
#include "tecplottools.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

char c1file[stringmaxlen];
char c2file[stringmaxlen];
char lofile[stringmaxlen];

int ntdens=10;
int ntpass=0;

class signed_curve *scurve1=NULL;
class signed_curve *scurve2=NULL;

int *I12=NULL;
int *I21=NULL;

int smart_dens=1;

double area(struct chain *l, double *w) {
  double a=(double)(0.0);
  double x,y,dx,dy;
  double dA;
  long i=0;
  struct chain *cur=l;
  while (cur->next!=NULL) {
    x=.5*(cur->x[0]+cur->next->x[0]);
    y=.5*(cur->x[1]+cur->next->x[1]);
    dx=cur->next->x[0]-cur->x[0];
    dy=cur->next->x[1]-cur->x[1];
    dA=y*dx-x*dy;
    if (w!=NULL) {
      dA*=w[i];
    }
    a+=0.5*dA;
    cur=cur->next;
    i++;
  }
  return(a);
}

double area(struct chain *l) {
  return(area(l,NULL));
}

void densify_local(int ndens, struct chain **c, long *n) {
  struct chain *cur=*c;
  struct chain *m;
  int i;
  double x,y;
  double x0,y0;
  double xf,yf;
  x0=cur->x[0];
  y0=cur->x[1];
  xf=cur->next->x[0];
  yf=cur->next->x[1];
  for (i=1;i<ndens;i++) {
    m=cur->next;
    cur->next=(struct chain *)malloc(sizeof(chain));
    cur->next->prev=cur;
    cur=cur->next;
    cur->next=m;
    cur->x[0]=x0+i*(xf-x0)/ndens;
    cur->x[1]=y0+i*(yf-y0)/ndens;
    *n=*n+1;
  }
  cur=cur->next;
  *c=cur;
}

void densCurves(int ndens, struct chain *c1, long *n, int sign2, struct chain *c2) {
  struct chain *cur=c1;
  int b1,b2;
  long i=0;
  signed_curve *scurve=NULL;
  scurve=chain2signedcurve(c2);
  I12=(int *)realloc(I12,*n*sizeof(int));
  fancy_display_reset("     Analyse  : ");
  while (cur!=NULL) {
    fancy_display_advance(i,*n);
    I12[i]=scurve->isInside(cur->x);
    i++;
    cur=cur->next;
  }
  fancy_display_stop();
  cur=c1;
  i=0;
  fancy_display_reset("     Densify  : ");
  while (cur->next!=NULL) {
    fancy_display_advance(i,*n);
    b1=I12[i];
    b2=I12[i+1];
    if ((b1==0)||(b2==0)) {
      i++;
      cur=cur->next;
    } else if (b1*b2<0) { 
      densify_local(ndens,&cur,n);
      i++;
    } else {
      i++;
      cur=cur->next;
    }
  }
  fancy_display_stop();
  delete scurve;
}

void densSmart(long *nc1, struct chain *c1, long *nc2, struct chain *c2) {
  struct chain *cur1=c1;
  struct chain *cur2=c2;
  struct chain *tmp;
  int b1,b2;
  long i=0;
  long j;
  double newX, newY;
  int stop;

  struct interOnCurve {
    struct chain *where;
    double p[2];
  };
  struct interOnCurve *inter=NULL;
  long ninter=0;

  signed_curve *s1=NULL;
  signed_curve *s2=NULL;
  s2=chain2signedcurve(c2);
  s1=chain2signedcurve(c1);
  I12=(int *)realloc(I12,*nc1*sizeof(int));
  I21=(int *)realloc(I21,*nc2*sizeof(int));
  fancy_display_reset("     Analyse A: ");
  while (cur1!=NULL) {
    fancy_display_advance(i,*nc1);
    I12[i]=s2->isInside(cur1->x);
    i++;
    cur1=cur1->next;
  }
  fancy_display_stop();
  cur1=c1;
  i=0;
  fancy_display_reset("     Analyse B: ");
  while (cur2!=NULL) {
    fancy_display_advance(i,*nc2);
    I21[i]=s1->isInside(cur2->x);
    i++;
    cur2=cur2->next;
  }
  fancy_display_stop();
  cur2=c2;
  i=0;
  fancy_display_reset("     Densifying:");
  while (cur1->next!=NULL) {
    fancy_display_advance(i,*nc1);
    b1=I12[i];
    b2=I12[i+1];
    if ((b1==0)||(b2==0)) {
      i++;
      cur1=cur1->next;
    } else if (b1*b2<0) { 
      //densify_local(ndens,&cur,n);
      cur2=c2;
      j=0;
      stop=0;
      while ((stop==0)&&(cur2->next!=NULL)) {
	if (I21[j]*I21[j+1]<=0) {
	  if (chainInter(cur1,cur2)) {
	    interCompute(cur1,cur2,&newX,&newY);
	    //printf("%lf\t%lf\n",newX,newY);
	    tmp=cur1->next;
	    cur1->next=(struct chain *)malloc(sizeof(struct chain));
	    cur1->next->prev=cur1;
	    cur1=cur1->next;
	    cur1->next=tmp;
	    cur1->x[0]=newX;
	    cur1->x[1]=newY;
	    tmp->prev=cur1;
	    *nc1=*nc1+1;
	    //remember the intersection point to densify B later:
	    ninter++;
	    inter=(struct interOnCurve *)realloc(inter,ninter*sizeof(struct interOnCurve));
	    inter[ninter-1].where=cur2;
	    inter[ninter-1].p[0]=newX;
	    inter[ninter-1].p[1]=newY;
	    stop=1;
	  } else {
	    j++;
	    cur2=cur2->next;
	  }
	} else {
	  j++;
	  cur2=cur2->next;
	}
      }
      i++;
      cur1=cur1->next;
    } else {
      i++;
      cur1=cur1->next;
    }
  }
  for (i=0;i<ninter;i++) {
    cur2=inter[i].where;
    tmp=cur2->next;
    cur2->next=(struct chain *)malloc(sizeof(struct chain));
    cur2->next->prev=cur2;
    cur2=cur2->next;
    cur2->next=tmp;
    cur2->x[0]=inter[i].p[0];
    cur2->x[1]=inter[i].p[1];
    tmp->prev=cur2;
    *nc2=*nc2+1;    
  }
  if (inter!=NULL) {
    free(inter);
  }
  fancy_display_stop();
  delete s1;
  delete s2;
  
}

/******
void densify(int ndens, struct chain *c, long *n) {
  struct chain *cur=c;
  struct chain *m;
  int i;
  double x,y;
  double x0,y0;
  double xf,yf;
  while (cur->next!=NULL) {
    x0=cur->x[0];
    y0=cur->x[1];
    xf=cur->next->x[0];
    yf=cur->next->x[1];
    for (i=1;i<ndens;i++) {
      m=cur->next;
      cur->next=(struct chain *)malloc(sizeof(chain));
      cur->next->prev=cur;
      cur=cur->next;
      cur->next=m;
      cur->x[0]=x0+i*(xf-x0)/ndens;
      cur->x[1]=y0+i*(yf-y0)/ndens;
      *n=*n+1;
    }
    cur=cur->next;
  }
}
********/


void centerg(struct chain *c, double *xg, double *yg) {
  struct chain *cur=c;
  double x,y,dx,dy,dA;
  double ar;
  *xg=0;
  *yg=0;
  while (cur->next!=NULL) {
    x=.5*(cur->x[0]+cur->next->x[0]);
    y=.5*(cur->x[1]+cur->next->x[1]);
    dx=cur->next->x[0]-cur->x[0];
    dy=cur->next->x[1]-cur->x[1];
    dA=-x*x*dy;
    *xg+=0.5*dA;
    dA=y*y*dx;
    *yg+=0.5*dA;
    cur=cur->next;
  }
  ar=area(c);
  *xg/=ar;
  *yg/=ar;
}

void enlarge(struct chain *c, double eps) {
  //first compute center of gravity
  double xg,yg;
  struct chain *cur;
  double v1,v2;
  double nv;
  //printf("TRANS/SHRINK/ENLARGE CURVE...\n");
  centerg(c,&xg,&yg);
  printf("  xg=%lf\tyg=%lf\n",xg,yg);
  cur=c;
  while (cur!=NULL) {
    v1=cur->x[0]-xg;
    v2=cur->x[1]-yg;
    nv=sqrt(v1*v1+v2*v2);
    cur->x[0]=cur->x[0]+eps*v1/nv;
    cur->x[1]=cur->x[1]+eps*v2/nv;
    cur=cur->next;
  }
}


int setWeights(int b1, int b2, double *q1, double *q2, double *q3) {

  switch (b1) {
  case 1:
    switch (b2) {
    case 1: *q1=-1.0; *q2=1.0; *q3=0.0; break;
    case 0: *q1=-1.0; *q2=1.0; *q3=0.0; break;
    case -1:*q1= 0.0; *q2=0.5; *q3=0.5; break;
    default:
      printf("LOBER-LIGHT-ERROR-FATAL\n");
      exit(1);
    }
    break;
  case 0:
    switch (b2) {
    case 1: *q1=-1.0; *q2=1.0; *q3=0.0; break;
    case 0: *q1= 0.0; *q2=0.5; *q3=0.5; break;
    case -1:*q1= 1.0; *q2=0.0; *q3=1.0; break;
    default:
      printf("LOBER-LIGHT-ERROR-FATAL\n");
      exit(1);
    }
    break;
  case -1:
    switch (b2) {
    case 1: *q1= 0.0; *q2=0.5; *q3=0.5; break;
    case 0: *q1= 1.0; *q2=0.0; *q3=1.0; break;
    case -1:*q1= 1.0; *q2=0.0; *q3=1.0; break;
    default:
      printf("LOBER-LIGHT-ERROR-FATAL\n");
      exit(1);
    }
    break;
  default:
    printf("LOBER-LIGHT-ERROR-FATAL\n");
    exit(1);
  }
  return(0);
}

void lober_light(int narg, char *arg[]) {
  struct chain *c1=NULL;
  struct chain *c2=NULL;
  struct chain *cur=NULL;
  long nc1, nc2;
  char *lline=NULL;
  double **w12=NULL;
  double **w21=NULL;
  long i,j;
  FILE *f;
  int b1,b2;

  FILE *f1,*f2,*f3,*f4, *f5, *f6;

  double arA, arB;
  int sarA,sarB;
  double arQ;
  double AinterB, AunionB;

  double gnrlerr;

  double in1,in2,out1,out2;
  double in,out,din,dout;

  strcpy(c1file,arg[2]);
  strcpy(c2file,arg[3]);
  strcpy(lofile,arg[4]);

  if (narg>=7) {
    smart_dens=0;
    sscanf(arg[6],"%i",&ntpass);
    ntdens=10;
    if (narg>=8) {
      sscanf(arg[7],"%i",&ntdens);
    }
  }

  printf("Loading the First curve from %s\n",c1file);
  f=fopen(c1file,"r");
  getLine(f,&lline);
  printf("  ***DataSetTitle=%sEOL\n",lline);
  nVar=TecTitle2nVar(lline);
  nc1=loadManiData(f,&c1);
  fclose(f);
  printf("%i points loaded\n",nc1);

  printf("Loading the Second curve from %s\n",c2file);
  f=fopen(c2file,"r");
  getLine(f,&lline);
  printf("  ***DataSetTitle=%sEOL\n",lline);
  nVar=TecTitle2nVar(lline);
  nc2=loadManiData(f,&c2);
  fclose(f);
  printf("%i points loaded\n",nc2);

  //enlarge(c1,1e-8);

  //printf("Translating to (xg,yg)...\n");
  //enlarge(c1,(double)(0.0));
  //enlarge(c2,(double)(0.0));

  if (smart_dens==1) {
    printf("Conditionning curves [1 smart pass]...\n");
    arA=area(c1);
    arB=area(c2);
    printf("[A]=%lf [B]=%lf nA=%i nB=%i\n",fabs(arA),fabs(arB),nc1,nc2);
    sarA=(arA>=(double)(0.0))?1:-1;
    sarB=(arB>=(double)(0.0))?1:-1;
    /***
    densCurves(1,c1,&nc1,sarB,c2);
    densCurves(1,c2,&nc2,sarA,c1);
    ***/
    //need to copy c1 here, so we can use the copy (no resampled) to
    //densify c2 afterwards !!!
    //TO DO
    densSmart(&nc1,c1,&nc2,c2);
    printf("[A]=%lf [B]=%lf nA=%i nB=%i\n",fabs(arA),fabs(arB),nc1,nc2);
    //densSmart(&nc2,c2,&nc1,c1);
    //printf("[A]=%lf [B]=%lf nA=%i nB=%i\n",fabs(arA),fabs(arB),nc1,nc2);
  } else {
    if ((ntpass>0)&&(ntdens>0)) {
      printf("Densifying curves [%ipass",ntpass);
      if (ntpass>1) {
	printf("es"); //What can I say? One step closer to perfection... 
      }
      printf("@%ipoint",ntdens);
      if (ntdens>1) {
	printf("s"); //What can I say? One step closer to perfection...
      }
      printf("]...\n",ntdens);
      for (i=0;i<ntpass;i++) {
	arA=area(c1);
	arB=area(c2);
	printf("[A]=%lf [B]=%lf nA=%i nB=%i\n",fabs(arA),fabs(arB),nc1,nc2);
	sarA=(arA>=(double)(0.0))?1:-1;
	sarB=(arB>=(double)(0.0))?1:-1;
	densCurves(ntdens,c1,&nc1,sarB,c2);
	densCurves(ntdens,c2,&nc2,sarA,c1);
      }
      printf("[A]=%lf [B]=%lf nA=%i nB=%i\n",fabs(arA),fabs(arB),nc1,nc2);
    }
  }

  printf("Checking data... ");
  arA=area(c1);
  arB=area(c2);
  sarA=(arA>=(double)(0.0))?1:-1;
  sarB=(arB>=(double)(0.0))?1:-1;
  if (arA*arB>=(double)(0.0)) {
    printf("OK!\n");
  } else {
    printf("FAIL!!!\n");
    printf("ERROR::Curves are not oriented in the same direction\n");
    printf("LOBER::Change direction on one curve or contact authors\n");
    printf("LOBER::EXITING NOW on ERROR\n");
    exit(1);
  }

  //here, I should create a "signed_curve" copy of each curve
  scurve1=chain2signedcurve(c1);
  scurve2=chain2signedcurve(c2);
  

  f1=fopen("c11.dat","w");
  fprintf(f1,"VARIABLES=\"x\"\"y\"\n");
  f2=fopen("c12.dat","w");
  fprintf(f2,"VARIABLES=\"x\"\"y\"\n");
  f3=fopen("c10.dat","w");
  fprintf(f3,"VARIABLES=\"x\"\"y\"\n");
  f4=fopen("c22.dat","w");
  fprintf(f4,"VARIABLES=\"x\"\"y\"\n");
  f5=fopen("c21.dat","w");
  fprintf(f5,"VARIABLES=\"x\"\"y\"\n");
  f6=fopen("c20.dat","w");
  fprintf(f6,"VARIABLES=\"x\"\"y\"\n");
  printf("Creating w_ij arrays...\n");
  w12=(double **)malloc(3*sizeof(double *));
  w21=(double **)malloc(3*sizeof(double *));
  for (i=0;i<3;i++) {
    w12[i]=(double *)malloc(nc1*sizeof(double));
    w21[i]=(double *)malloc(nc2*sizeof(double));
  }
  I12=(int *)realloc(I12,nc1*sizeof(int));
  I21=(int *)realloc(I21,nc2*sizeof(int));
  cur=c1;
  i=0;
  I12[0]=scurve2->isInside(cur->x);
  fancy_display_reset("     W12 :      ");
  while (cur->next!=NULL) {
    fancy_display_advance(i,nc1);
    //b1=isInside_wrapper_chain(sarB,cur->x[0],cur->x[1],c2);
    //b2=isInside_wrapper_chain(sarB,cur->next->x[0],cur->next->x[1],c2);
    //b1=scurve2->isInside(cur->x);
    b1=I12[i];
    b2=scurve2->isInside(cur->next->x);
    I12[i+1]=b2;
    if (b1==1) {
      fprintf(f1,OUT_FORMAT,cur->x[0]);
      fprintf(f1,"\t");
      fprintf(f1,OUT_FORMAT,cur->x[1]);
      fprintf(f1,"\n");
    } else if (b1==-1) {
      fprintf(f2,OUT_FORMAT,cur->x[0]);
      fprintf(f2,"\t");
      fprintf(f2,OUT_FORMAT,cur->x[1]);
      fprintf(f2,"\n");
    } else {
      fprintf(f3,OUT_FORMAT,cur->x[0]);
      fprintf(f3,"\t");
      fprintf(f3,OUT_FORMAT,cur->x[1]);
      fprintf(f3,"\n");
    }

    setWeights(b1,b2,&(w12[0][i]),&(w12[1][i]),&(w12[2][i]));

    /*******
    if ((b1==1)&&(b2==1)) {
      w12[0][i]=(double)(-1.0);
      w12[1][i]=(double)(1.0);
      w12[2][i]=(double)(0.0);
    } else if ((b1==-1)&&(b2==-1)) {
      w12[0][i]=(double)(1.0);
      w12[1][i]=(double)(0.0);
      w12[2][i]=(double)(1.0);
    }
    } else {
      w12[0][i]=(double)(0.0);
      w12[1][i]=(double)(0.5);
      w12[2][i]=(double)(0.5);
    }
    ********/
    i++;
    cur=cur->next;
  }
  fancy_display_stop();
  cur=c2;
  i=0;
  I21[0]=scurve1->isInside(cur->x);
  fancy_display_reset("     W21 :      ");
  while (cur->next!=NULL) {
    fancy_display_advance(i,nc2);
    //b1=isInside_wrapper_chain(sarA,cur->x[0],cur->x[1],c1);
    //b2=isInside_wrapper_chain(sarA,cur->next->x[0],cur->next->x[1],c1);
    //b1=scurve1->isInside(cur->x);
    b1=I21[i];
    b2=scurve1->isInside(cur->next->x);
    I21[i+1]=b2;
    if (b1==1) {
      fprintf(f4,OUT_FORMAT,cur->x[0]);
      fprintf(f4,"\t");
      fprintf(f4,OUT_FORMAT,cur->x[1]);
      fprintf(f4,"\n");
    } else if (b1==-1) {
      fprintf(f5,OUT_FORMAT,cur->x[0]);
      fprintf(f5,"\t");
      fprintf(f5,OUT_FORMAT,cur->x[1]);
      fprintf(f5,"\n");
    } else {
      fprintf(f6,OUT_FORMAT,cur->x[0]);
      fprintf(f6,"\t");
      fprintf(f6,OUT_FORMAT,cur->x[1]);
      fprintf(f6,"\n");
    }

    setWeights(b1,b2,&(w21[0][i]),&(w21[1][i]),&(w21[2][i]));

    /**********
    if (b1&&b2) {
      w21[0][i]=(double)(-1.0);
      w21[1][i]=(double)(1.0);
      w21[2][i]=(double)(0.0);
    } else if ((b1==0)&&(b2==0)) {
      w21[0][i]=(double)(1.0);
      w21[1][i]=(double)(0.0);
      w21[2][i]=(double)(1.0);
    } else {
      w21[0][i]=(double)(0.0);
      w21[1][i]=(double)(0.5);
      w21[2][i]=(double)(0.5);
    }
    ***********/
    i++;
    cur=cur->next;
  }
  fancy_display_stop();
  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);
  fclose(f5);
  fclose(f6);

  arQ=sarA*area(c1,w12[0]);
  arQ+=sarB*area(c2,w21[0]);
  AinterB=sarA*area(c1,w12[1]);
  AinterB+=sarB*area(c2,w21[1]);
  AunionB=sarA*area(c1,w12[2]);
  AunionB+=sarB*area(c2,w21[2]);

  arA=fabs(arA);
  arB=fabs(arB);

  printf("[A]=  %lf\n",arA);
  printf("[B]=  %lf\n",arB);
  printf("[Q]=  %lf\n",arQ);
  printf("[AXB]=%lf\n",AinterB);
  printf("[AUB]=%lf\n",AunionB);

  AunionB=fabs(AunionB);
  AinterB=fabs(AinterB);
  arQ=fabs(arQ);

  in=.25*(AunionB-AinterB+arQ)+.5*(arA-arB);
  out=.25*(AunionB-AinterB+arQ)+.5*(arB-arA);
  din=.5*fabs(AunionB-AinterB-arQ);
  dout=din;
  din/=in;
  dout/=out;

  gnrlerr=(AunionB+AinterB-arA-arB);
  gnrlerr=fabs(gnrlerr);
  gnrlerr/=(arA+arB);

  printf("***RESULTS SECTION**********\n");
  printf("   LOBE_IN= %lf (error=%e%%)\n",in,100*din);
  printf("   LOBE_OUT=%lf (error=%e%%)\n",out,100*dout);
  printf("   AREA_ERROR=%e%%\n",100*gnrlerr);

  f=fopen(lofile,"w");
  fprintf(f,OUT_FORMAT,in);
  fprintf(f,"\t");
  fprintf(f,OUT_FORMAT,100*din);
  fprintf(f,"\t");
  fprintf(f,OUT_FORMAT,out);
  fprintf(f,"\t");
  fprintf(f,OUT_FORMAT,100*dout);
  fprintf(f,"\n");
  fclose(f);

  for (i=0;i<3;i++) {
    free(w12[i]);
    free(w21[i]);
  }
  free(w12);
  free(w21);
  if (I12!=NULL) {
    free(I12);
  }
  if (I21!=NULL) {
    free(I21);
  }

  delete scurve1;
  delete scurve2;

}

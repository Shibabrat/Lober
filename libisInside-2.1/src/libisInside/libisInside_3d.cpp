#undef  __ISINSIDE__3D__
#define __ISINSIDE__3D__INTG
#undef  __ISINSIDE__3D__INTGINTG

#include <isInside.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void signed_surface::make_finalize(void) {
  int k,l;
  this->majorPt=(int *)malloc(ntri*sizeof(int));
  this->majorV1xV2=(double *)malloc(ntri*sizeof(double));
  for (k=0;k<2;k++) {
    this->V1[k]=(double *)malloc(ntri*sizeof(double));
    this->V2[k]=(double *)malloc(ntri*sizeof(double));
  }
  for (k=0;k<3;k++) {
    this->normal[k]=(double *)malloc(ntri*sizeof(double));
    for (l=0;l<3;l++) {
      this->rotation[k][l]=(double *)malloc(ntri*sizeof(double));
    }
  }
  this->myVolume=volume_compute();
  this->rot=1;
  if (this->myVolume<(double)(0.0)) {
    this->myVolume*=-1;
    this->rot=-1;
  }
}

void signed_surface::make(long N, long M,double *x, long *tri) {
  long i;
  int k;
  this->npt=N;
  this->ntri=M;
  for (k=0;k<3;k++) {
    this->pt[k]=(double *)malloc(npt*sizeof(double));
    this->tri[k]=(long *)malloc(ntri*sizeof(long));
  }
  for (i=0;i<N;i++) {
    for (k=0;k<3;k++) {
      this->pt[k][i]=x[3*i+k];
    }
  }
  for (i=0;i<M;i++) {
    for (k=0;k<3;k++) {
      this->tri[k][i]=tri[3*i+k];
    }
  }
  make_finalize();
}

signed_surface::signed_surface(long N, long M,double *x, long *tri) {
  this->make(N,M,x,tri);
}

signed_surface::~signed_surface() {
  int k,l;
  free(majorPt);
  free(majorV1xV2);
  for (k=0;k<2;k++) {
    free(V1[k]);
    free(V2[k]);
  }
  for (k=0;k<3;k++) {
    free(pt[k]);
    free(tri[k]);
    free(normal[k]);
    for (l=0;l<3;l++) {
      free(rotation[k][l]);
    }
  }
}

double signed_surface::volume(void) {
  return(myVolume);
}

void signed_surface::normal_vector_tryout(long faceID, int majorPoint, double *v1, double *v2, double *v1xv2_norm, double v1xv2[3], double ve1[3], double ve2[3]) {
  //double ve1[3];
  //double ve2[3];
  int p1, p2, p3;
  int k;
  //for a face #faceID, this function starts at point #majorPoint
  // I.E. 0<=majorPoint<=2. It construct the vector v1 and v2 from majorPoint
  // to the other vertex of the face in such a way that v1xv2 is oriented
  // in the positive direction of the normal vector to the face.
  //The function returns ||v1||, ||v2|| and ||v1xv2||.
  //This function is used to try majorPoint=0,1,2 and pick which one will
  //give the best numerical result for the normal vector (v1xv2)/||v1xv2||
  p1=majorPoint;
  p2=p1+1;
  p2=(p2>2)?0:p2;
  p3=p2+1;
  p3=(p3>2)?0:p3;
  for (k=0;k<3;k++) {
    ve1[k]=pt[k][tri[p2][faceID]]-pt[k][tri[p1][faceID]];
    ve2[k]=pt[k][tri[p3][faceID]]-pt[k][tri[p1][faceID]];
  }
  v1xv2[0]=ve1[1]*ve2[2]-ve2[1]*ve1[2];
  v1xv2[1]=ve1[2]*ve2[0]-ve2[2]*ve1[0];
  v1xv2[2]=ve1[0]*ve2[1]-ve2[0]*ve1[1];
  //*v1v2=(double)(0.0);
  *v1=(double)(0.0);
  *v2=(double)(0.0);
  *v1xv2_norm=(double)(0.0);
  for (k=0;k<3;k++) {
    *v1+=ve1[k]*ve1[k];
    *v2+=ve2[k]*ve2[k];
    *v1xv2_norm+=v1xv2[k]*v1xv2[k];
    //*v1v2+=ve1[k]*ve2[k];
  }
  *v1=sqrt(*v1);
  *v2=sqrt(*v2);
  *v1xv2_norm=sqrt(*v1xv2_norm);
  //printf("n=%lf\t%lf\t%lf\n",v1xv2[0],v1xv2[1],v1xv2[2]);
}

double signed_surface::volume_compute_face(long i) {
  double V=(double)(0.0);
  double v1,v2,v1xv2_n,v1xv2[3];
  double b_v1,b_v2,b_v1xv2_n,b_v1xv2[3];
  int biangle=0;
  double b_angle;
  double angle;
  double teta, alpha, proj;
  double ve1[3], ve2[3], b_ve1[3], b_ve2[3];
  int k,l;
  //double tmp1[3], tmp2[3];
  this->normal_vector_tryout(i,biangle,&b_v1,&b_v2,&b_v1xv2_n,b_v1xv2,b_ve1,b_ve2);
  if ((b_v1<=(double)(0.0))||(b_v2<=(double)(0.0))) {
    //this face has no surface !
    return((double)(0.0));
  }
  b_angle=fabs(b_v1xv2_n/(b_v1*b_v2));
  for (k=1;k<3;k++) {
    this->normal_vector_tryout(i,k,&v1,&v2,&v1xv2_n,v1xv2,ve1,ve2);
    if ((v1<=(double)(0.0))||(v2<=(double)(0.0))) {
      //this face has no surface !
      return((double)(0.0));
    }
    angle=fabs(v1xv2_n/(v1*v2));
    if (angle>b_angle) {
      biangle=k;
      b_v1=v1;
      b_v2=v2;
      b_v1xv2_n=v1xv2_n;
      //b_v1v2=v1v2;
      for (l=0;l<3;l++) {
	b_v1xv2[l]=v1xv2[l];
	b_ve1[l]=ve1[l];
	b_ve2[l]=ve2[l];
      }
      b_angle=angle;
    }
  }
  //printf("%i:%i\n",i,biangle);
  //save some of this data for future use
  majorPt[i]=biangle;
  majorV1xV2[i]=b_v1xv2_n;
  for (k=0;k<3;k++) {
    V+=(pt[k][tri[biangle][i]]*b_v1xv2[k])/((double)(6.0));
  }
  //here we set the normal and rotation data for this face
  for (k=0;k<3;k++) {
    normal[k][i]=b_v1xv2[k];
  }
  teta=atan2(normal[1][i],normal[0][i]);
  proj=sqrt(normal[0][i]*normal[0][i]+normal[1][i]*normal[1][i]);
  if (proj<=(double)(0.0)) {
    alpha=M_PI/((double)(2.0));
    if (normal[2][i]>=(double)(0.0)) {
      alpha*=-1;
    }
  } else {
    alpha=-atan(normal[2][i]/proj);
  }
  //printf("%i:%lf\t%lf\t%lf\n",i,normal[0][i],normal[1][i],normal[2][i]);
  //printf("%i:%lf\t%lf\t%lf\n",i,b_ve1[0],b_ve1[1],b_ve1[2]);
  //printf("%i:%lf\t%lf\t%lf\n",i,b_ve2[0],b_ve2[1],b_ve2[2]);
  //printf("%i:teta=%lf\talpha=%lf\n",i,teta,alpha);
  //printf("V=%lf\n",V);
  rotation[0][0][i]=cos(teta)*cos(alpha);
  rotation[0][1][i]=sin(teta)*cos(alpha);
  rotation[0][2][i]=-sin(alpha);
  rotation[1][0][i]=-sin(teta);
  rotation[1][1][i]=cos(teta);
  rotation[1][2][i]=(double)(0.0);
  rotation[2][0][i]=cos(teta)*sin(alpha);
  rotation[2][1][i]=sin(teta)*sin(alpha);
  rotation[2][2][i]=cos(alpha);

  //Rotate ve1 and ve2 and save data for future use
  for (k=1;k<3;k++) {
    //tmp1[k]=(double)(0.0);
    //tmp2[k]=(double)(0.0);
    V1[k-1][i]=(double)(0.0);
    V2[k-1][i]=(double)(0.0);
    for (l=0;l<3;l++) {
      V1[k-1][i]+=rotation[k][l][i]*b_ve1[l];
      V2[k-1][i]+=rotation[k][l][i]*b_ve2[l];
    }
  }
  //printf("v1:%lf\t%lf\n",V1[0][i],V1[1][i]);
  //printf("v2:%lf\t%lf\n",V2[0][i],V2[1][i]);

  /***
  printf("%i:",i);
  double t;
  for (k=0;k<3;k++) {
    t=(double)(0.0);
    for (l=0;l<3;l++) {
      t+=rotation[k][l][i]*normal[l][i];
    }
    printf("%lf\t",t);
  }
  printf("\n");
  ***/
  return(V);
}

double signed_surface::volume_compute(void) {
  //this function computes the (signed) volume contained in the surface
  double V=(double)(0.0);
  long i;
  for (i=0;i<ntri;i++) {
    //first we search for the pair of vectors on the face that
    //have an angle closest to 90 degrees (best normal vector)
    V+=volume_compute_face(i);
  }
  return(V);
}

double signed_surface::d_isInside(double x[3]) {
  double ret=(double)(0.0);
  long i;
  int bndFlag=0;
  for (i=0;i<ntri;i++) {
    ret+=this->d_isInside_face(i,x,&bndFlag);
    if (bndFlag!=0) {
      return(2*M_PI);
    }
  }
  return(rot*ret);
}

#ifdef __ISINSIDE__3D__INTGINTG
double signed_surface::d_isInside_face(long i, double x[3], int *bndFlag) {
  double ret=(double)(0.0);
  double dret;
  double xrot[3];
  int k, l;
  double r;
  //double r0=(double)(0.0);
  long nt, ns;
  long nT=10;
  long nS=10;
  double ds, dt;
  double dy, dz;
  *bndFlag=0; //not on the boundary
  for (k=0;k<3;k++) {
    xrot[k]=(double)(0.0);
    for (l=0;l<3;l++) {
      xrot[k]+=rotation[k][l][i]*(pt[l][tri[majorPt[i]][i]]-x[l]);
    }
    //r0+=xrot[k]*xrot[k];
  }
  
  ds=1.0/nS;
  dt=1.0/nT;
  for (nt=0;nt<nT;nt++) {
    for (ns=0;ns<(1-nt*dt)/ds;ns++) {
      dy=V1[0][i]*nt*dt+V2[0][i]*ns*ds;
      dz=V1[1][i]*nt*dt+V2[1][i]*ns*ds;
      r=xrot[0]*xrot[0]+(xrot[1]+dy)*(xrot[1]+dy)+(xrot[2]+dz)*(xrot[2]+dz);
      r=sqrt(r);
      if (r<=(double)(0.0)) {
	*bndFlag=1; //majorPoint on the boundary
	return((double)(0.0));
      }
      dret=xrot[0]*majorV1xV2[i]/r;
      dret/=r;
      dret/=r;
      ret+=ds*dt*dret;
    }
  }
  return(ret);
}
#endif

#ifdef __ISINSIDE__3D__INTG
double signed_surface::d_isInside_face(long i, double x[3], int *bndFlag) {
  double ret=(double)(0.0);
  double dret;
  double r;
  long nt, ns;
  long nT=10;
  long nS=10;
  double dt;
  double dy, dz;

  int k, l;
  double P,Q,R,delta,eta,Delta;
  double xrot[3];
  double t;

  *bndFlag=0; //not on the boundary
  for (k=0;k<3;k++) {
    xrot[k]=(double)(0.0);
    for (l=0;l<3;l++) {
      xrot[k]+=rotation[k][l][i]*(pt[l][tri[majorPt[i]][i]]-x[l]);
    }
  }

  P=V1[0][i]*V2[0][i]+V1[1][i]*V2[1][i];
  Q=V1[0][i]*xrot[1]+V1[1][i]*xrot[2];
  R=V2[0][i]*xrot[1]+V2[1][i]*xrot[2];
  delta=V1[0][i]*V1[0][i]+V1[1][i]*V1[1][i];
  eta=V2[0][i]*V2[0][i]+V2[1][i]*V2[1][i];
  Delta=xrot[0]*xrot[0]+xrot[1]*xrot[1]+xrot[2]*xrot[2];
  delta=sqrt(delta);
  eta=sqrt(eta);
  Delta=sqrt(Delta);

  dt=1.0/nT;

  double p1, p2, p3;

  for (nt=0;nt<nT;nt++) {
    t=nt*dt;
    p1=t*(P-eta*eta)+R+eta*eta;
    p2=t*t*(delta*delta*eta*eta-P*P)+2*t*(Q*eta*eta-R*P)+eta*eta*Delta*Delta-R*R;
    p3=t*t*(delta*delta-eta*eta-2*P)+2*t*(P+Q-R-eta*eta)+eta*eta+2*R+Delta*Delta;
    dret=sqrt(-p3);



    ret+=dret;
  }

  /***
  for (nt=0;nt<nT;nt++) {
    for (ns=0;ns<(1-nt*dt)/ds;ns++) {
      dy=V1[0][i]*nt*dt+V2[0][i]*ns*ds;
      dz=V1[1][i]*nt*dt+V2[1][i]*ns*ds;
      r=xrot[0]*xrot[0]+(xrot[1]+dy)*(xrot[1]+dy)+(xrot[2]+dz)*(xrot[2]+dz);
      r=sqrt(r);
      if (r<=(double)(0.0)) {
	*bndFlag=1; 
	return((double)(0.0));
      }
      dret=xrot[0]*majorV1xV2[i]/r;
      dret/=r;
      dret/=r;
      ret+=ds*dt*dret;
    }
  }
  ****/

  return(ret);
}
#endif

#include <isInside.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_ARR 100

long tri[MAX_ARR][3];
double pt[MAX_ARR][3];
long npt=0;
long ntri=0;

#define xmin -1.5
#define xmax  1.5
#define nx    30

void export_mesh(const char *fname) {
  long i;
  FILE *f;
  f=fopen(fname,"w");
  if (f==NULL) {
    system("mkdir data");
    f=fopen(fname,"w");
  }
  if (f==NULL) {
    fprintf(stderr,"Cannot open %s\n",fname);
    exit(1);
  }

  fprintf(f,"VARIABLES=\"x\"\"y\"\"z\"\n");
  fprintf(f,"ZONE N=%i, E=%i, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",npt,ntri);
  for (i=0;i<npt;i++) {
    fprintf(f,"%lf\t%lf\t%lf\n",pt[i][0],pt[i][1],pt[i][2]);
  }
  for (i=0;i<ntri;i++) {
    fprintf(f,"%i\t%i\t%i\n",tri[i][0]+1,tri[i][1]+1,tri[i][2]+1);
  }
  fclose(f);
}

double shape(double tet, double phi) {
  return(1.0);
}

void reshape_mesh(void) {
  long i;
  double x,y,z;
  double tet, phi;
  double r;
  for (i=0;i<npt;i++) {
    x=pt[i][0];
    y=pt[i][1];
    z=pt[i][2];
    r=sqrt(x*x+y*y+z*z);
    phi=atan2(y,x);
    tet=asin(z/r);
    r=shape(tet,phi);
    pt[i][0]=r*cos(tet)*cos(phi);
    pt[i][1]=r*cos(tet)*sin(phi);
    pt[i][2]=r*sin(tet);
  }
}

void divide_faces(void) {
  long rem_nfaces=ntri;
  long i;
  int k;
  long p1, p2, p3;
  for (i=0;i<rem_nfaces;i++) {
    //divide face i
    npt++;
    for (k=0;k<3;k++) {
      pt[npt-1][k]=(pt[tri[i][0]][k]+pt[tri[i][1]][k]+pt[tri[i][2]][k])/3.0;
    }
    p1=tri[i][0];p2=tri[i][1];p3=tri[i][2];
    tri[i][2]=npt-1;
    ntri++;
    tri[ntri-1][0]=p2;    tri[ntri-1][1]=p3;    tri[ntri-1][2]=npt-1;
    ntri++;
    tri[ntri-1][0]=p3;    tri[ntri-1][1]=p1;    tri[ntri-1][2]=npt-1;
  }
}

int main(int narg, char **arg) {

  long i;
  int k;
  signed_surface *s;
  double arr_pt[3*MAX_ARR];
  long arr_tri[3*MAX_ARR];
  long a,b,c;
  double x[3];
  FILE *f;

  npt=8;
  pt[0][0]=-1; pt[0][1]=-1; pt[0][2]=-1;
  pt[1][0]= 1; pt[1][1]=-1; pt[1][2]=-1;
  pt[2][0]=-1; pt[2][1]= 1; pt[2][2]=-1;
  pt[3][0]= 1; pt[3][1]= 1; pt[3][2]=-1;
  pt[4][0]=-1; pt[4][1]=-1; pt[4][2]= 1;
  pt[5][0]= 1; pt[5][1]=-1; pt[5][2]= 1;
  pt[6][0]=-1; pt[6][1]= 1; pt[6][2]= 1;
  pt[7][0]= 1; pt[7][1]= 1; pt[7][2]= 1;
  ntri=12;
  tri[0][0]=0;  tri[0][1]=1;  tri[0][2]=2;
  tri[1][0]=2;  tri[1][1]=1;  tri[1][2]=3;
  tri[2][0]=2;  tri[2][1]=3;  tri[2][2]=6;
  tri[3][0]=3;  tri[3][1]=7;  tri[3][2]=6;
  tri[4][0]=6;  tri[4][1]=7;  tri[4][2]=4;
  tri[5][0]=7;  tri[5][1]=5;  tri[5][2]=4;
  tri[6][0]=3;  tri[6][1]=1;  tri[6][2]=5;
  tri[7][0]=5;  tri[7][1]=7;  tri[7][2]=3;
  tri[8][0]=1;  tri[8][1]=0;  tri[8][2]=4;
  tri[9][0]=4;  tri[9][1]=5;  tri[9][2]=1;
  tri[10][0]=0; tri[10][1]=2; tri[10][2]=6;
  tri[11][0]=6; tri[11][1]=4; tri[11][2]=0;
  
  //reshape_mesh();
  //divide_faces();
  //reshape_mesh();

  export_mesh("data/boundary3D.dat");

  for (i=0;i<npt;i++) {
    for (k=0;k<3;k++) {
      arr_pt[3*i+k]=pt[i][k];
    }
  }
  for (i=0;i<ntri;i++) {
    for (k=0;k<3;k++) {
      arr_tri[3*i+k]=tri[i][k];
    }
  }
  s=new signed_surface(npt,ntri,arr_pt,arr_tri);

  printf("Volume=%lf\n",s->volume());

  f=fopen("data/test3D.dat","w");
  if (f==NULL) {
    system("mkdir data");
    f=fopen("data/test3D.dat","w");
  }
  if (f==NULL) {
    fprintf(stderr,"Cannot open %s\n","data/test3D.dat");
    exit(1);
  }

  fprintf(f,"VARIABLES=\"x\"\"y\"\"z\"\"d\"\n");
  fprintf(f,"ZONE I=%i, J=%i, K=%i\n",nx+1,nx+1,nx+1);

  for (a=0;a<=nx;a++) {
    x[0]=xmin+a*(xmax-xmin)/nx;
    for (b=0;b<=nx;b++) {
      x[1]=xmin+b*(xmax-xmin)/nx;
      for (c=0;c<=nx;c++) {
	x[2]=xmin+c*(xmax-xmin)/nx;
	fprintf(f,"%lf\t%lf\t%lf\t%lf\n",x[0],x[1],x[2],s->d_isInside(x));
      }
    }
  }

  fclose(f);

  delete s;

  return(0);
}

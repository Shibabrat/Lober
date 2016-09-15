/*********************************************/
/* LOBER 1.5                                 */
/*    (C) Francois Lekien, Shane Ross        */
/*        California Inst. of Technology     */
/*********************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <isInside.h>
#include "lober.h"
#include "tecplottools.h"
#include "inandout.h"
#include "envVar.h"
#include "lober_light.h"
#include "compatibility.h"
#include "config.h"


char unsfile[stringmaxlen];
char stbfile[stringmaxlen];
char manifile[stringmaxlen];
char lobefile[stringmaxlen];
char pipsfile[stringmaxlen];
char lobebfile[stringmaxlen];


/***DATASET PARAMS***
#define unsskip 300
#define stbskip 100
#define zones  500
/***END***/

/***TEST***
#define unsskip 500
#define stbskip 0
#define zones  601
/***END***/

int unsskip=-1;
int stbskip=-1;
int zones=-1;
int lobeseq=-1;

/***TEST***
#define unsskip 0
#define stbskip 0
#define zones  1
/***END***/


int nVar=0;
int nx, ny;
int light_mode=0;
int nbip=-1;

void help_line() {
  printf("usage1: lober <unsfile> <stbfile> <lobefile>\n");
  printf("              <manifile> <pipsfile> unsskip\n");
  printf("               stbskip nZones [<lobebfile>[lobeseq]]\n");
  printf("               [-BIP nbib]\n");
  printf("description: unsfile = file containing unstable manifold\n");
  printf("             stbfile = file containing stable manifold\n");
  printf("             lobefile = file to output colored lobes\n");
  printf("             manifile = file to output simplified manifolds\n");
  printf("             pipsfile = file to output pips\n");
  printf("             unsskip, stbskip = zones to skip in uns/stbfile\n");
  printf("             nZones = number of zones\n");
  printf("             lobebfile = file to output lobe boundary\n");
  printf("             lobeseq = divide lobes in sequences (default=2)\n");
  printf("             -BIP nbip = which pip defines the boundary\n");
  printf("                         (starting at \"0\" ; default=npips/2)\n");
  printf("usage2: lober -light <separfile> <lobefile> <lobeoutfile>\n");
  printf("              [ -DENS nPass nDens ]\n");
  printf("description: separfile = file containing separatrix\n");
  printf("             lobefile = files containing lobes\n");
  printf("             lobeoutfile = files to output cleaned lobes\n");
  printf("             -DENS : densify curves based on intersections\n");
  printf("                     nPass = number of passes\n");
  printf("                     nDens = number of points to add\n");
  exit(0);
}

void arg_error(const char *txt) {
  printf("lober::invalid command line\n");
  printf("lober::%s\n",txt);
  printf("lober::try lober -h for help\n");
  exit(1);
}

void parseArg(int narg, char *arg[]) {
  int i;
  int shift=0;
  for (i=0;i<narg;i++) {
    if ((strcmp(arg[i],"-h")==0)||(strcmp(arg[i],"--help")==0)) {
      help_line();
    }
    if (strcmp(arg[i],"-light")==0) {
      light_mode=1;
    }
    if ((nbip<0)&&(strcmp(arg[i],"-BIP")==0)) {
      nbip=i;
    }
  }

  if (nbip>=0) {
    if (nbip<narg-1) {
      if (sscanf(arg[nbip+1],"%i",&nbip)<=0) {
	arg_error("-BIP must be followed by positive integer");
      }
      if (nbip<0) {
	arg_error("-BIP must be followed by positive integer");
      }
    } else {
      arg_error("-BIP must be followed by positive integer");
    }
  }

  shift=(nbip<0)?0:2;

  if (light_mode==1) {
    if ((narg!=5)&&(narg!=8)&&(narg!=7)) {
      arg_error("invalid number of arguments in light mode");
    }
    if ((narg==8)||(narg==7)) {
      if (strcmp(arg[5],"-DENS")!=0) {
	arg_error("syntax error in light mode");
      }
    }
    return;
  }

  if ((narg<9+shift)||(narg>11+shift)) {
    arg_error("invalid number of arguments");
  }

  strcpy(unsfile,arg[1]);
  strcpy(stbfile,arg[2]);
  strcpy(lobefile,arg[3]);
  strcpy(manifile,arg[4]);
  strcpy(pipsfile,arg[5]);
  if (sscanf(arg[6],"%i",&unsskip)!=1)
    arg_error("unsskip must be an integer!");
  if (sscanf(arg[7],"%i",&stbskip)!=1)
    arg_error("stbskip must be an integer!");
  if (sscanf(arg[8],"%i",&zones)!=1)
    arg_error("nZones must be an integer!");

  if (narg>9+shift) {
    strcpy(lobebfile,arg[9]);
    lobeseq=2; //by default
    if (narg>10+shift) {
      if (sscanf(arg[10],"%i",&lobeseq)!=1) {
	arg_error("lobeseq must be an integer!");
      }
      if (lobeseq<1) {
	arg_error("lobeseq must be strictly positive!\n");
      }
    }
  }

  if ((unsskip<0)||(stbskip<0))
    arg_error("unsskip and stbskip must be >=0!");

  if (zones<=0)
    arg_error("nZones must be > 0!");
  
}

int stat_is_inside(void) {
  unsigned int v1, v2, v3;
  
  printf("LOBER: loading isInside library...");
  fflush(stdout);
  v1=__IS_INSIDE__VERSION_INTERFACE__();
  v2=__IS_INSIDE__VERSION_REVISION__();
  v3=__IS_INSIDE__VERSION_COMPATIBILITY__();
  printf(" lib%s v%i.%i.%i\n",__IS_INSIDE__NAME__(),v1,v2,v3);
  if (v1<COMPAT_IS_INSIDE_INTERFACE) {
    printf("       LOBER requires libisInside 2.x or above\n");
    printf("LOBER: Please upgrade libisInside to its latest version!\n");
    printf("       email: lekien AT princeton DOT edu\n");
    printf("       web: http://www.lekien.com/~francois/software/libisInside\n");
    printf("LOBER: Abnormal termination!\n");
    exit(1);
  }
  if (v1-v3>COMPAT_IS_INSIDE_INTERFACE) {
    printf("       ERROR: libisInside %i.%i.%i is not compatible with its 2.x interface\n",v1,v2,v3);
    printf("LOBER: The version of libisInside installed on this system is\n");
    printf("       not compatible with libisInside 2.x that was used to\n");
    printf("       build this version of Lober.\n");
    printf("LOBER: Please upgrade Lober to its latest version!\n");
    printf("       email: lekien AT princeton DOT edu\n");
    printf("       web: http://www.lekien.com/~francois/software/lober\n");
    printf("LOBER: Abnormal termination!\n");
    exit(1);
  }
}

int main(int narg, char *arg[]) {

	long i;
	int j,k;
	double vx,vy;
	long *stbpos;
	FILE *stbin, *unsin;
	struct chain *stb=NULL;
	struct chain *uns=NULL;
	long nstb, nuns;
	FILE *maniin, *lobein;
	FILE *pipsin;
	FILE *lobebin;
	struct pip *pips=NULL;
	int npips;
	int nsep;
	int nlobes=0;
	struct lobe **lobes=NULL;
	struct lobe *separ=NULL;
	struct pip *firstpip;
	struct pip *bpip;
	struct pip *lastpip;
	struct pip *p;
	char omode[2048];
	signed_curve **lobe_curves=NULL;

	printf("LOBER: lober v1.8 03-15-2005\n");
	printf("       Francois Lekien ( lekien AT princeton DOT edu )\n");
	printf("       Shane Ross ( shane AT cds DOT caltech DOT edu )\n");
	
	printf("           running on %s (%s)\n",getStringOSName(),getStringCName());
	printf("           user=%s group=%s\n",getUserName(),getUserGroup());
	printf("           pwd=%s\n",getPath());

	parseArg(narg, arg);

	stat_is_inside();

	if (light_mode==1) {
	  lober_light(narg,arg);
	} else {

	stbpos=(long *)malloc(zones*sizeof(long));

	printf("marking stable manifold...\n");
	markstbmanifold(zones,stbpos,stbfile);
	printf("stable manifold marked!\n");

	printf("loading unstable manifold file...\n");
	unsin=setunsmanifold(unsfile);
	printf("unstable manifold loaded\n");
	stbin=fopen(stbfile,"r");
	maniin=fopen(manifile,"w");
	lobein=fopen(lobefile,"w");
	pipsin=fopen(pipsfile,"w");

	if (maniin==NULL) {
	  printf("LOBER::Unable to open %s\n",manifile);
	  exit(0);
	}
	if (lobein==NULL) {
	  printf("LOBER::Unable to open %s\n",lobefile);
	  exit(0);
	}
	if (pipsin==NULL) {
	  printf("LOBER::Unable to open %s\n",pipsfile);
	  exit(0);
	}

	if (lobeseq>0) {
	  lobebin=fopen(lobebfile,"w");
	  if (lobebin==NULL) {
	    printf("LOBER::Unable to open %s\n",lobebfile);
	    exit(0);
	  }
	} else {
	  lobebin=NULL;
	}

	fprintf(maniin,"VARIABLES=\"x\"\"y\"\n");
	fprintf(lobein,"VARIABLES=\"x\"\"y\"\"lobe\"\n");
	fprintf(pipsin,"VARIABLES=\"x\"\"y\"\n");
	if (lobebin!=NULL) {
	  fprintf(lobebin,"VARIABLES=\"x\"\"y\"\n");
	}
	
	printf("entering the main loop\n");

	for (i=0;i<zones;i++) {
		eraseChain(stb);
		eraseChain(uns);
		stb=NULL;
		uns=NULL;
		nuns=loadManiData(unsin,&uns);
		//add a first point
		//addMani(-10.0,uns->x[1],&uns);
		//nuns++;
		fseek(stbin,stbpos[i],SEEK_SET);
		nstb=loadManiData(stbin,&stb);
		//addMani(2010.0,stb->x[1],&stb);
		//nstb++;

		/***here we get the lobes***/

		pips=getPips(uns,stb,&npips);
		printf("There are %i pips\n",npips);

		lobes=getLobes(npips,pips,&nlobes);
		printf("There are %i lobes\n",nlobes);

		/***output***/
		printf("output...\n");
		firstpip=pips;
		lastpip=pips;
		if (pips!=NULL) {
			while (lastpip->next!=NULL) {
				lastpip=lastpip->next;
			}
		}
		if ((cutEnd)||(firstpip==NULL)) {
		  //nstb=decimateMani(decimdmin,stb,firstpip->stb[1]);
		  //nuns=decimateMani(decimdmin,uns,lastpip->uns[1]);
		} else {
		  //nstb=decimateMani(decimdmin,stb,NULL);
		  //nuns=decimateMani(decimdmin,uns,NULL);
		}
		fprintf(maniin,"ZONE I=%i\n",nuns);
		outputMani(maniin,uns);
		fprintf(maniin,"ZONE I=%i\n",nstb);
		outputMani(maniin,stb);
		fprintf(pipsin,"ZONE I=%i\n",npips);
		outputPips(pipsin,pips);

		/***DEBUGCODE
		//for (k=0;k<nlobes;k++) {
		k=8;
			outputLobe(lobein,lobes[k]);
		//}
		exit(0);
		***ENDDEBUG***/

		//this is where I should decimate the lobes!
		for (j=0;j<nlobes;j++) {
		  //decimateLobe(decimdmin,lobes[j]);
		}


		if (lobeseq>0) {
		  outputLobeb(lobebin,i,nlobes,lobes);
		}

		sprintf(omode,"%s\t%s\t%%i\n",OUT_FORMAT,OUT_FORMAT);

		lobe_curves=(signed_curve **)malloc(nlobes*sizeof(signed_curve *));
		for (k=0;k<nlobes;k++) {
		  lobe_curves[k]=lobe2signedcurve(lobes[k]);
		}
		fprintf(lobein,"ZONE I=%i J=%i\n",lbnx,lbny);
		for (k=0;k<lbny;k++) {
			vy=ymin+k*(ymax-ymin)/(lbny-1);
			for (j=0;j<lbnx;j++) {
				vx=xmin+j*(xmax-xmin)/(lbnx-1);
				//fprintf(lobein,omode,vx,vy,isInLobes(vx,vy,nlobes,lobes));
				fprintf(lobein,omode,vx,vy,isInLobes(vx,vy,nlobes,lobe_curves));
			}
		}
		for (k=0;k<nlobes;k++) {
		  delete lobe_curves[k];
		}
		free(lobe_curves);

		printf("npips=%i\n",npips);

		/*********************************************/
		/* This is where we define the BIP           */
		/*********************************************/
		nsep=(npips-1)/2;
		if (nbip>=0) {
		  if (nbip<npips) {
		    nsep=nbip;
		    printf("Using user input for nbip=%i\n",nbip);
		  } else {
		    printf("lober::Error User specified nbip=%i\n",nbip);
		    printf("     ::But there are only %i pips\n",npips);
		    printf("lober::REVISE THE ARG of -BIP nbip and re-run\n");
		    exit(0);
		  }
		}
		printf("Building separatrix at pip=%i\n",nsep);
		bpip=firstpip;
		for (j=0;j<nsep;j++) {
		  bpip=bpip->next;
		}
		/*********************************************/
		/* End of BIP Definition                     */
		/*********************************************/
		separ=getSepar(uns,stb,bpip);
		outputSepar(separ);
		getLobeStat(separ);
		
		p=firstpip;
		printf("LOBE\tAREA\tIPs\n");
		for (j=0;j<nlobes;j++) {
		  printf("%3i ",j+1);
		  getLobeStat(lobes[j],p,p->next,separ,bpip);
		  p=p->next;
		}

		printf("freeing...\n");
		if (pips!=NULL)
			erasePips(pips);

		for (j=0;j<nlobes;j++) {
			if (lobes[j]!=NULL) {
				eraseLobe(lobes[j]);
			}
		}
		if (lobes!=NULL) {
			free(lobes);
		}

	}

	fclose(stbin);
	fclose(unsin);
	fclose(lobein);
	fclose(maniin);
	if (lobebin!=NULL) {
	  fclose(lobebin);
	}
	free(stbpos);
	}

	printf("LOBER::DONE!\n");
	
	return(0);
}

void outputSepar(struct lobe *sep) {
  FILE *f;
  struct lobe *c=sep;
  char omode[2048];
  f=fopen("separ.dat","w");
  fprintf(f,"VARIABLES=\"x\"\"y\"\n");
  sprintf(omode,"%s\t%s\n",OUT_FORMAT,OUT_FORMAT);
  while (c!=NULL) {
    fprintf(f,omode,c->x[0],c->x[1]);
    c=c->next;
  }
  fclose(f);
}

struct lobe *getSepar(struct chain *u, struct chain *s, struct pip *bip) {
  struct lobe *sep=NULL;
  struct lobe *ret=NULL;
  struct chain *c;
  c=u;
  while (c!=bip->uns[1]) {
    if (sep==NULL) {
      sep=(struct lobe *)malloc(sizeof(struct lobe));
      ret=sep;
    } else {
      sep->next=(struct lobe *)malloc(sizeof(struct lobe));
      sep=sep->next;
    }
    sep->x[0]=c->x[0];
    sep->x[1]=c->x[1];
    sep->next=NULL;
    sep->rot=1;
    c=c->next;
  }
  //printf("h1\n");
  sep->next=(struct lobe *)malloc(sizeof(struct lobe));
  sep=sep->next;
  sep->x[0]=bip->x[0];
  sep->x[1]=bip->x[1];
  sep->next=NULL;
  sep->rot=1;
  c=bip->stb[0];
  //printf("h2\n");
  while (c!=s) {
    sep->next=(struct lobe *)malloc(sizeof(struct lobe));
    sep=sep->next;
    sep->x[0]=c->x[0];
    sep->x[1]=c->x[1];
    sep->next=NULL;
    sep->rot=1;
    c=c->prev;
  }
  //printf("h3\n");
  sep->next=(struct lobe *)malloc(sizeof(struct lobe));
  sep=sep->next;
  sep->x[0]=c->x[0];
  sep->x[1]=c->x[1];
  sep->next=NULL;
  sep->rot=1;
  //printf("h4\n");
  return(ret);
}

void addMani(double x, double y, struct chain **c) {
	//add an extra point at the beginning
	struct chain *p;
	p=*c;
	*c=(struct chain *)malloc(sizeof(struct chain));
	(*c)->x[0]=x;
	(*c)->x[1]=y;
	(*c)->prev=NULL;
	(*c)->next=p;
	p->prev=(*c);
}

void decimateLobe(double dmin, struct lobe *l) {
	long npts=0;
	long del=0;
	int needDel;
	struct lobe *cur=l;
	struct lobe *rm;
	double d;

	while (cur->next!=NULL) {
		npts++;
		needDel=0;
		d=(cur->x[0]-cur->next->x[0])*(cur->x[0]-cur->next->x[0])+(cur->x[1]-cur->next->x[1])*(cur->x[1]-cur->next->x[1]);
		if (d<dmin*dmin) {
			needDel=1;
		}
		if (needDel) {
			if (cur->next->next!=NULL) {
				//delete point
				del++;
				rm=cur->next;
				cur->next=rm->next;
				free(rm);
			} else {
				cur=cur->next;
			}
		} else {
			cur=cur->next;
		}
	}
	
	printf("decimator:%i->%i\n",npts,npts-del);
}

long decimateMani(double dmin, struct chain *l, struct chain *theEnd) {
	long npts=0;
	long del=0;
	int needDel;
	struct chain *cur=l;
	struct chain *rm;
	double d;

	while (cur->next!=theEnd) {
		npts++;
		needDel=0;
		d=(cur->x[0]-cur->next->x[0])*(cur->x[0]-cur->next->x[0])+(cur->x[1]-cur->next->x[1])*(cur->x[1]-cur->next->x[1]);
		if (d<dmin*dmin) {
			needDel=1;
		}
		if (needDel) {
			if (cur->next->next!=NULL) {
				//delete point
				del++;
				rm=cur->next;
				cur->next=rm->next;
				free(rm);
			} else {
				cur=cur->next;
			}
		} else {
			cur=cur->next;
		}
	}

	
	cur=cur->next; //cur=theEnd;
	npts++;

	if (cur!=NULL) {
	while (cur->next!=NULL) {
		npts++;
		del++;
		rm=cur->next;
		cur->next=rm->next;
		free(rm);	
	}
	npts++;
	}

//	npts++;

	printf("decimatorMANI:%i->%i\n",npts,npts-del);
	return(npts-del);
}

void outputLobeb(FILE *f, long zonen, int nlobes, struct lobe **l) {
  int i;
  struct lobe *cur;
  char *fname;
  char tmp[100];
  int j;
  long ireal=0;
  long npoints;
  char omode[2048];
  FILE *fil2;
  fname=(char *)malloc((strlen(lobebfile)+100)*sizeof(char));
  for (i=0;i<nlobes;i++) {
    npoints=0;
    cur=l[i];
    while (cur!=NULL) {
      cur=cur->next;
      npoints++;
    }
    if (npoints>=LOBEMINNUMPOINTS) {
      strcpy(fname,lobebfile);
      sprintf(tmp,".%5i",ireal);
      for (j=0;j<strlen(tmp);j++) {
	if (tmp[j]==' ') {
	  tmp[j]='0';
	}
      }
      strcat(fname,tmp);
      //printf("%s\n",fname);
      fil2=fopen(fname,"w");
      fprintf(f,"ZONE T=\"zone=%i lobe=%i\"\n",zonen,ireal);
      cur=l[i];
      sprintf(omode,"%s\t%s\n",OUT_FORMAT,OUT_FORMAT);
      while (cur!=NULL) {
	fprintf(f,omode,cur->x[0],cur->x[1]);
	fprintf(fil2,omode,cur->x[0],cur->x[1]);
	cur=cur->next;
      }
      fclose(fil2);
      ireal++;
    }
  }
  if (ireal!=nlobes) {
    printf("There were only %i lobes printed\n",ireal);
  }
  free(fname);
}

int isInLobes(double cx, double cy, int nlobes, class signed_curve **l) {
  int val=0;
  int i=0;
  int isI;
  while (i<nlobes) {
    isI=isInLobe(cx,cy,l[i]);
    if (isI!=0) {
      val=isI;
      i=nlobes;
    }
    i++;
  }
  return(val);
}

int isInLobe(double cx, double cy, class signed_curve *l) {
  int ret=0;
  int isI;
  double p[2];
  p[0]=cx;
  p[1]=cy;
  //return(isInside_wrapper_lobe(cx,cy,l));
  isI=l->isInside(p);
  if ((isI==1)||(isI==0)) {
    ret=l->rotation();
  }
  return(ret);
}

struct lobe **getLobes(int npips, struct pip *pips, int *nlobes) {
	struct lobe **l=NULL;
	struct pip *cur;
	int i;
	*nlobes=npips-1;
	if (*nlobes<=0) {
		return(NULL);
	}
	l=(struct lobe **)malloc(*nlobes*sizeof(struct lobe *));
	cur=pips;
	for (i=0;i<*nlobes;i++) {
		//printf("Lobe %i\n",i+1);
		l[i]=getLobe(cur);
		cur=cur->next;

	}
	return(l);
}

double area(struct lobe *l) {
  double a=(double)(0.0);
  double x,y,dx,dy;
  double dA;
  struct lobe *cur=l;
  while (cur->next!=NULL) {
    x=.5*(cur->x[0]+cur->next->x[0]);
    y=.5*(cur->x[1]+cur->next->x[1]);
    dx=cur->next->x[0]-cur->x[0];
    dy=cur->next->x[1]-cur->x[1];
    dA=y*dx-x*dy;
    a+=0.5*dA;
    cur=cur->next;
  }
  return(a);
}

void getLobeStat(struct lobe *l) {
  getLobeStat(l,NULL,NULL,NULL,NULL);
}

void getLobeStat(struct lobe *l, struct pip *p1, struct pip *p2, struct lobe *sep, struct pip *bip) {
  long nsips1=0;
  long nsips2=0;
  long i;
  long n1=0,n2=0;
  struct chain *c;
  struct chain *s;
  struct pip **ip1=NULL;
  struct pip **ip2=NULL;
  struct lobe *subl=NULL;
  struct lobe *subli=NULL;
  struct chain *start, *stop;
  int mode;
  int sw,fnd;
  double a_in=(double)(0.0);
  double a_out=(double)(0.0);
  char lmode[2048];

  sprintf(lmode,"\t%s",STDOUT_FORMAT);
  printf(lmode,fabs(area(l)));
  if ((p1!=NULL)&&(p2!=NULL)&&(sep!=NULL)&&(bip!=NULL)) {

    c=p1->uns[0];
    while (c!=p2->uns[1]) {
      s=bip->stb[1];
      while (s->prev!=NULL) {
	if (chainInter(c,s->prev)) {
	  nsips1++;
	  ip1=(struct pip **)realloc(ip1,nsips1*sizeof(struct pip *));
	  ip1[nsips1-1]=(struct pip *)malloc(sizeof(struct pip));
	  ip1[nsips1-1]->uns[0]=c;
	  ip1[nsips1-1]->uns[1]=c->next;
	  ip1[nsips1-1]->stb[0]=s->prev;
	  ip1[nsips1-1]->stb[1]=s;
	  ip1[nsips1-1]->next=NULL;
	  if (nsips1>1) {
	    ip1[nsips1-1]->prev=ip1[nsips1-2];
	    ip1[nsips1-2]->next=ip1[nsips1-1];
	  } else {
	    ip1[nsips1-1]->prev=NULL;
	  }
	  ip1[nsips1-1]->rot=1;
	  interCompute(c,s->prev,&(ip1[nsips1-1]->x[0]),&(ip1[nsips1-1]->x[1]));	  
	}
	s=s->prev;
      }
      c=c->next;
    }
    printf("\t%i",nsips1);

    c=p2->stb[0];
    while (c!=p1->stb[1]) {
      s=bip->uns[1];
      while (s->prev!=NULL) {
	if (chainInter(c,s->prev)) {
	  nsips2++;
	  ip2=(struct pip **)realloc(ip2,nsips2*sizeof(struct pip *));
	  ip2[nsips2-1]=(struct pip *)malloc(sizeof(struct pip));
	  ip2[nsips2-1]->uns[0]=c;
	  ip2[nsips2-1]->uns[1]=c->next;
	  ip2[nsips2-1]->stb[0]=s->prev;
	  ip2[nsips2-1]->stb[1]=s;
	  ip2[nsips2-1]->next=NULL;
	  if (nsips2>1) {
	    ip2[nsips2-1]->prev=ip2[nsips2-2];
	    ip2[nsips2-2]->next=ip2[nsips2-1];
	  } else {
	    ip2[nsips2-1]->prev=NULL;
	  }
	  ip2[nsips2-1]->rot=1;
	  interCompute(c,s->prev,&(ip2[nsips2-1]->x[0]),&(ip2[nsips2-1]->x[1]));	  
	}
	s=s->prev;
      }
      c=c->next;
    }

    printf(",%i",nsips2);
  

    mode=0;
    if (nsips1>1) {
      mode+=1;
    }
    if (nsips2>1) {
      mode+=2;
    }
    
    if (((nsips1!=1)&&(nsips1%2==1))||((nsips2!=1)&&(nsips2%2==1))) {
      mode=5;
      printf("LOBE GEOM. MADNESS ");
    }
    
    if ((mode<1)||(mode>2)) {
      printf("::nsips1=%i nsips2=%i",nsips1,nsips2);
    }
    
    //add the first pip at the end
    
    if (mode==2) {
      //sub area computation mode==2
      nsips2++;
      ip2=(struct pip **)realloc(ip2,nsips2*sizeof(struct pip *));
      ip2[nsips2-1]=(struct pip *)malloc(sizeof(struct pip));
      ip2[nsips2-1]->uns[0]=ip2[0]->uns[0];
      ip2[nsips2-1]->uns[1]=ip2[0]->uns[1];
      ip2[nsips2-1]->stb[0]=ip2[0]->stb[0];
      ip2[nsips2-1]->stb[1]=ip2[0]->stb[1];
      ip2[nsips2-1]->next=NULL;
      if (nsips2>1) {
	ip2[nsips2-1]->prev=ip2[nsips2-2];
	ip2[nsips2-2]->next=ip2[nsips2-1];
      } else {
	ip2[nsips2-1]->prev=NULL;
      }
      ip2[nsips2-1]->rot=ip2[0]->rot;
      ip2[nsips2-1]->x[0]=ip2[0]->x[0];
      ip2[nsips2-1]->x[1]=ip2[0]->x[1];
      
      sw=0;
      
      for (i=0;i<nsips2-1;i++) {
	//printf("i=%i\n",i);
	if (sw==0) {
	  //integrate along stable manifold to next point
	  //we don't need to determine the direction first
	  //but we do it to crash from numercal errors
	  //start=ip2[i]->uns[1];
	  //stop=ip2[i+1]->uns[1];
	  start=ip2[i]->stb[1];
	  stop=ip2[i+1]->stb[1];
	  fnd=0;
	  while ((start!=NULL)&&(start!=stop)) {
	    start=start->next;
	  }
	  if (start==NULL) {
	    fnd=-1;
	  } else {
	    fnd=1;
	  }
	  start=ip2[i]->stb[1];

	  if (subl==NULL) {
	    subl=(struct lobe *)malloc(sizeof(struct lobe));
	    subli=subl;
	  } else {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	  }
	  subl->next=NULL;
	  subl->x[0]=ip2[i]->x[0];
	  subl->x[1]=ip2[i]->x[1];
	  subl->rot=1;
	  while (start!=stop) {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	    subl->next=NULL;
	    subl->x[0]=start->x[0];
	    subl->x[1]=start->x[1];
	    subl->rot=1;
	    if (fnd>0) {
	      start=start->next;
	    } else {
	      start=start->prev;
	    }
	  }
	  subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	  subl=subl->next;
	  subl->next=NULL;
	  subl->x[0]=ip2[i+1]->x[0];
	  subl->x[1]=ip2[i+1]->x[1];
	  subl->rot=1;
	  sw=1;
	} else {
	  //integrate along the unstable manifold
	  //we need to determine the direction first
	  start=ip2[i]->uns[1];
	  stop=ip2[i+1]->uns[1];
	  fnd=0;
	  while ((start!=NULL)&&(start!=stop)) {
	    start=start->next;
	  }
	  if (start==NULL) {
	    fnd=-1;
	  } else {
	    fnd=1;
	  }
	  start=ip2[i]->uns[1];
	  //now the loop
	  subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	  subl=subl->next;
	  subl->next=NULL;
	  subl->x[0]=ip2[i]->x[0];
	  subl->x[1]=ip2[i]->x[1];
	  subl->rot=1;
	  while (start!=stop) {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	    subl->next=NULL;
	    subl->x[0]=start->x[0];
	    subl->x[1]=start->x[1];
	    subl->rot=1;
	    if (fnd<0) {
	      start=start->prev;
	    } else {
	      start=start->next;
	    }
	  }
	  subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	  subl=subl->next;
	  subl->next=NULL;
	  subl->x[0]=ip2[i+1]->x[0];
	  subl->x[1]=ip2[i+1]->x[1];
	  subl->rot=1;
	  sw=0;
	}
	//printf("done\n");
      }
      a_in=area(subli);
    }

    if (mode==1) {
      //sub area computation mode==1
      nsips1++;
      ip1=(struct pip **)realloc(ip1,nsips1*sizeof(struct pip *));
      ip1[nsips1-1]=(struct pip *)malloc(sizeof(struct pip));
      ip1[nsips1-1]->uns[0]=ip1[0]->uns[0];
      ip1[nsips1-1]->uns[1]=ip1[0]->uns[1];
      ip1[nsips1-1]->stb[0]=ip1[0]->stb[0];
      ip1[nsips1-1]->stb[1]=ip1[0]->stb[1];
      ip1[nsips1-1]->next=NULL;
      if (nsips1>1) {
	ip1[nsips1-1]->prev=ip1[nsips1-2];
	ip1[nsips1-2]->next=ip1[nsips1-1];
      } else {
	ip1[nsips1-1]->prev=NULL;
      }
      ip1[nsips1-1]->rot=ip1[0]->rot;
      ip1[nsips1-1]->x[0]=ip1[0]->x[0];
      ip1[nsips1-1]->x[1]=ip1[0]->x[1];
      
      sw=0;
      
      for (i=0;i<nsips1-1;i++) {
	//printf("i=%i\n",i);
	if (sw==0) {
	  //integrate along unstable manifold to next point
	  //we don't need to determine the direction first
	  //but we do it to crash from numercal errors
	  //start=ip2[i]->uns[1];
	  //stop=ip2[i+1]->uns[1];
	  start=ip1[i]->uns[1];
	  stop=ip1[i+1]->uns[1];
	  fnd=0;
	  while ((start!=NULL)&&(start!=stop)) {
	    start=start->next;
	  }
	  if (start==NULL) {
	    fnd=-1;
	  } else {
	    fnd=1;
	  }
	  start=ip1[i]->uns[1];

	  if (subl==NULL) {
	    subl=(struct lobe *)malloc(sizeof(struct lobe));
	    subli=subl;
	  } else {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	  }
	  subl->next=NULL;
	  subl->x[0]=ip1[i]->x[0];
	  subl->x[1]=ip1[i]->x[1];
	  subl->rot=1;
	  while (start!=stop) {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	    subl->next=NULL;
	    subl->x[0]=start->x[0];
	    subl->x[1]=start->x[1];
	    subl->rot=1;
	    if (fnd>0) {
	      start=start->next;
	    } else {
	      start=start->prev;
	    }
	  }
	  subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	  subl=subl->next;
	  subl->next=NULL;
	  subl->x[0]=ip1[i+1]->x[0];
	  subl->x[1]=ip1[i+1]->x[1];
	  subl->rot=1;
	  sw=1;
	} else {
	  //integrate along the stable manifold
	  //we need to determine the direction first
	  start=ip1[i]->stb[1];
	  stop=ip1[i+1]->stb[1];
	  fnd=0;
	  while ((start!=NULL)&&(start!=stop)) {
	    start=start->next;
	  }
	  if (start==NULL) {
	    fnd=-1;
	  } else {
	    fnd=1;
	  }
	  start=ip1[i]->stb[1];
	  //now the loop
	  subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	  subl=subl->next;
	  subl->next=NULL;
	  subl->x[0]=ip1[i]->x[0];
	  subl->x[1]=ip1[i]->x[1];
	  subl->rot=1;
	  while (start!=stop) {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	    subl->next=NULL;
	    subl->x[0]=start->x[0];
	    subl->x[1]=start->x[1];
	    subl->rot=1;
	    if (fnd<0) {
	      start=start->prev;
	    } else {
	      start=start->next;
	    }
	  }
	  subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	  subl=subl->next;
	  subl->next=NULL;
	  subl->x[0]=ip1[i+1]->x[0];
	  subl->x[1]=ip1[i+1]->x[1];
	  subl->rot=1;
	  sw=0;
	}
	//printf("done\n");
      }
      a_out=area(subli);
    }

    if ((mode>=1)&&(mode<=2)) {
      
      subl=subli;
      while (subl!=NULL) {
	subli=subl;
	subl=subl->next;
	free(subli);
      }
      subl=NULL;
      subli=NULL;
    }


    //now compute the other area

    if (mode==2) {
      //sub area computation mode==2

      sw=1;
      
      for (i=0;i<nsips2-1;i++) {
	//printf("i=%i\n",i);

	if (i==nsips2-2) {
	  //if this is the last sip, stay on unstable manifold!
	  sw=1;
	}

	if (sw==0) {
	  //integrate along stable manifold to next point
	  //we don't need to determine the direction first
	  //but we do it to crash from numercal errors
	  //start=ip2[i]->uns[1];
	  //stop=ip2[i+1]->uns[1];
	  start=ip2[i]->stb[1];
	  stop=ip2[i+1]->stb[1];
	  fnd=0;
	  while ((start!=NULL)&&(start!=stop)) {
	    start=start->next;
	  }
	  if (start==NULL) {
	    fnd=-1;
	  } else {
	    fnd=1;
	  }
	  start=ip2[i]->stb[1];

	  if (subl==NULL) {
	    subl=(struct lobe *)malloc(sizeof(struct lobe));
	    subli=subl;
	  } else {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	  }
	  subl->next=NULL;
	  subl->x[0]=ip2[i]->x[0];
	  subl->x[1]=ip2[i]->x[1];
	  subl->rot=1;
	  while (start!=stop) {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	    subl->next=NULL;
	    subl->x[0]=start->x[0];
	    subl->x[1]=start->x[1];
	    subl->rot=1;
	    if (fnd>0) {
	      start=start->next;
	    } else {
	      start=start->prev;
	    }
	  }
	  subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	  subl=subl->next;
	  subl->next=NULL;
	  subl->x[0]=ip2[i+1]->x[0];
	  subl->x[1]=ip2[i+1]->x[1];
	  subl->rot=1;
	  sw=1;
	} else {
	  //integrate along the unstable manifold
	  //we need to determine the direction first
	  start=ip2[i]->uns[1];
	  stop=ip2[i+1]->uns[1];
	  fnd=0;
	  while ((start!=NULL)&&(start!=stop)) {
	    start=start->next;
	  }
	  if (start==NULL) {
	    fnd=-1;
	  } else {
	    fnd=1;
	  }
	  start=ip2[i]->uns[1];
	  //now the loop

	  if (subl==NULL) {
	    subl=(struct lobe *)malloc(sizeof(struct lobe));
	    subli=subl;
	  } else {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	  }

	  //subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	  //subl=subl->next;
	  subl->next=NULL;
	  subl->x[0]=ip2[i]->x[0];
	  subl->x[1]=ip2[i]->x[1];
	  subl->rot=1;
	  while (start!=stop) {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	    subl->next=NULL;
	    subl->x[0]=start->x[0];
	    subl->x[1]=start->x[1];
	    subl->rot=1;
	    if (fnd<0) {
	      start=start->prev;
	    } else {
	      start=start->next;
	    }
	  }
	  subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	  subl=subl->next;
	  subl->next=NULL;
	  subl->x[0]=ip2[i+1]->x[0];
	  subl->x[1]=ip2[i+1]->x[1];
	  subl->rot=1;
	  sw=0;
	}
	//printf("done\n");
      }
      a_out=area(subli);
    }

    if (mode==1) {
      //sub area computation mode==1
            
      sw=1;
      
      for (i=0;i<nsips1-1;i++) {
	//printf("i=%i\n",i);

	if (i==nsips1-2) {
	  //if this is the last sip, stay on stable manifold!
	  sw=1;
	}

	if (sw==0) {
	  //integrate along unstable manifold to next point
	  //we don't need to determine the direction first
	  //but we do it to crash from numercal errors
	  //start=ip2[i]->uns[1];
	  //stop=ip2[i+1]->uns[1];
	  start=ip1[i]->uns[1];
	  stop=ip1[i+1]->uns[1];
	  fnd=0;
	  while ((start!=NULL)&&(start!=stop)) {
	    start=start->next;
	  }
	  if (start==NULL) {
	    fnd=-1;
	  } else {
	    fnd=1;
	  }
	  start=ip1[i]->uns[1];

	  if (subl==NULL) {
	    subl=(struct lobe *)malloc(sizeof(struct lobe));
	    subli=subl;
	  } else {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	  }
	  subl->next=NULL;
	  subl->x[0]=ip1[i]->x[0];
	  subl->x[1]=ip1[i]->x[1];
	  subl->rot=1;
	  while (start!=stop) {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	    subl->next=NULL;
	    subl->x[0]=start->x[0];
	    subl->x[1]=start->x[1];
	    subl->rot=1;
	    if (fnd>0) {
	      start=start->next;
	    } else {
	      start=start->prev;
	    }
	  }
	  subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	  subl=subl->next;
	  subl->next=NULL;
	  subl->x[0]=ip1[i+1]->x[0];
	  subl->x[1]=ip1[i+1]->x[1];
	  subl->rot=1;
	  sw=1;
	} else {
	  //integrate along the stable manifold
	  //we need to determine the direction first
	  start=ip1[i]->stb[1];
	  stop=ip1[i+1]->stb[1];
	  fnd=0;
	  while ((start!=NULL)&&(start!=stop)) {
	    start=start->next;
	  }
	  if (start==NULL) {
	    fnd=-1;
	  } else {
	    fnd=1;
	  }
	  start=ip1[i]->stb[1];
	  //now the loop

	  if (subl==NULL) {
	    subl=(struct lobe *)malloc(sizeof(struct lobe));
	    subli=subl;
	  } else {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	  }
	  //subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	  //subl=subl->next;
	  subl->next=NULL;
	  subl->x[0]=ip1[i]->x[0];
	  subl->x[1]=ip1[i]->x[1];
	  subl->rot=1;
	  while (start!=stop) {
	    subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	    subl=subl->next;
	    subl->next=NULL;
	    subl->x[0]=start->x[0];
	    subl->x[1]=start->x[1];
	    subl->rot=1;
	    if (fnd<0) {
	      start=start->prev;
	    } else {
	      start=start->next;
	    }
	  }
	  subl->next=(struct lobe *)malloc(sizeof(struct lobe));
	  subl=subl->next;
	  subl->next=NULL;
	  subl->x[0]=ip1[i+1]->x[0];
	  subl->x[1]=ip1[i+1]->x[1];
	  subl->rot=1;
	  sw=0;
	}
	//printf("done\n");
      }
      a_in=area(subli);
    }

    if ((mode>=1)&&(mode<=2)) {
      
      subl=subli;
      while (subl!=NULL) {
	subli=subl;
	subl=subl->next;
	free(subli);
      }
    }

    sprintf(lmode,"\t%s\t%s\t%s",STDOUT_FORMAT,STDOUT_FORMAT,STDOUT_FORMAT);
    if ((mode>=1)&&(mode<=2)) {
      printf(lmode,fabs(a_in),fabs(a_out),fabs(a_in+a_out));
    }

  }
  
  for (i=0;i<nsips1;i++) {
    free(ip1[i]);
  }
  if (ip1!=NULL) {
    free(ip1);
  }
  for (i=0;i<nsips2;i++) {
    free(ip2[i]);
  }
  if (ip2!=NULL) {
    free(ip2);
  }
  
  printf("\n");
  
}

void getLobeStat_old(struct lobe *l, struct pip *p1, struct pip *p2, struct lobe *sep) {
  long nsips=0;
  struct lobe *c;
  struct lobe *s;
  struct chain *vc1;
  struct chain *vc2;
  struct chain *vc3;
  struct chain *vc4;
  struct lobe **ip_c=NULL;
  struct lobe **ip_s=NULL;
  struct lobe *tmp1=NULL;
  struct lobe *tmp2=NULL;
  char lmode[2048];

  vc1=(struct chain *)malloc(sizeof(struct chain));
  vc2=(struct chain *)malloc(sizeof(struct chain));
  vc1->next=vc2;
  vc2->prev=vc1;
  vc1->prev=NULL;
  vc2->next=NULL;
  vc3=(struct chain *)malloc(sizeof(struct chain));
  vc4=(struct chain *)malloc(sizeof(struct chain));
  vc3->next=vc4;
  vc4->prev=vc3;
  vc3->prev=NULL;
  vc4->next=NULL;
  sprintf(lmode,"LOBE area=%s\n",STDOUT_FORMAT);
  printf(lmode,area(l));
  if ((p1!=NULL)&&(p2!=NULL)&&(sep!=NULL)) {
    nsips=1;
    ip_c=(struct lobe **)realloc(ip_c,sizeof(struct lobe *));
    ip_s=(struct lobe **)realloc(ip_s,sizeof(struct lobe *));
    ip_c[nsips-1]=l;
    ip_s[nsips-1]=NULL;
    c=l->next;
    while (c->next->next!=NULL) {
      vc3->x[0]=c->x[0];
      vc3->x[1]=c->x[1];
      vc4->x[0]=c->next->x[0];
      vc4->x[1]=c->next->x[1];

      s=sep;
      while (s->next!=NULL) {
	vc1->x[0]=s->x[0];
	vc1->x[1]=s->x[1];
	vc2->x[0]=s->next->x[0];
	vc2->x[1]=s->next->x[1];
	
	
	if (chainInter(vc3,vc1)) {
	  nsips++;
	  //we need to add the ip
	  tmp1=c->next;
	  tmp2=s->next;
	  c->next=(struct lobe *)malloc(sizeof(struct lobe));
	  s->next=(struct lobe *)malloc(sizeof(struct lobe));
	  interCompute(vc3,vc1,&(c->next->x[0]),&(c->next->x[1]));
	  s->next->x[0]=c->next->x[0];
	  s->next->x[1]=c->next->x[1];
	  s->next->rot=s->rot;
	  c->next->rot=c->rot;

	  ip_c=(struct lobe **)realloc(ip_c,nsips*sizeof(struct lobe *));
	  ip_s=(struct lobe **)realloc(ip_s,nsips*sizeof(struct lobe *));
	  ip_c[nsips-1]=c->next;
	  ip_s[nsips-1]=s->next;

	  s->next->next=tmp2;
	  c->next->next=tmp1;
	  s=s->next;
	  c=c->next;
	}
	
	s=s->next;
      }
      c=c->next;
    }

    nsips++;
    ip_c=(struct lobe **)realloc(ip_c,nsips*sizeof(struct lobe *));
    ip_s=(struct lobe **)realloc(ip_s,nsips*sizeof(struct lobe *));
    ip_c[nsips-1]=ip_c[0];
    ip_s[nsips-1]=ip_s[0];

    printf("nsips=%i\n",nsips);
  }
  free(vc1);
  free(vc2);
  free(vc3);
  free(vc4);
  if (ip_c!=NULL) {
    free(ip_c);
  }
  if (ip_s!=NULL) {
    free(ip_s);
  }
}

struct lobe *getLobe(struct pip *p) {
	struct lobe *l;
	struct lobe *cur;
	struct chain *c;
	double rottmp;

	//add the first IP
	l=(struct lobe *)malloc(sizeof(struct lobe));
	l->x[0]=p->x[0];
	l->x[1]=p->x[1];
	l->next=NULL;
	l->rot=p->rot;
	//printf("pip1:%f %f\n",p->x[0],p->x[1]);
	//printf("pip2:%f %f\n",p->next->x[0],p->next->x[1]);
	//add the segment of uns
	cur=l;
	c=p->uns[1];
	
	//if (c!=p->next->uns[1]) {
	  while (c!=p->next->uns[1]) {
	    
	    cur->next=(struct lobe *)malloc(sizeof(struct lobe));
	    cur=cur->next;
	    cur->x[0]=c->x[0];
	    cur->x[1]=c->x[1];
	    cur->next=NULL;
	    cur->rot=p->rot;
	    
	    c=c->next;
	  }
	  //}
	
	//add the second IP
	cur->next=(struct lobe *)malloc(sizeof(struct lobe));
	cur=cur->next;
	cur->x[0]=p->next->x[0];
	cur->x[1]=p->next->x[1];
	//printf("pip3:%f %f\n",cur->x[0],cur->x[1]);
	//printf("pip4:%f %f\n",p->next->x[0],p->next->x[1]);
	cur->next=NULL;
	cur->rot=p->rot;
	//add the segment of stb
	/*
	c=p->stb[1];
	if (c!=p->next->stb[1]) {
	while (c!=p->next->stb[0]) {
	*/	
	c=p->next->stb[1];
	//if (c!=p->stb[0]) {
	while (c!=p->stb[1]) {

		cur->next=(struct lobe *)malloc(sizeof(struct lobe));
		cur=cur->next;
		cur->x[0]=c->x[0];
		cur->x[1]=c->x[1];
		cur->next=NULL;
		cur->rot=p->rot;
		
		c=c->next;
	}
	//}
	//add the first IP
	cur->next=(struct lobe *)malloc(sizeof(struct lobe));
	cur=cur->next;
	cur->x[0]=p->x[0];
	cur->x[1]=p->x[1];
	cur->next=NULL;
	cur->rot=p->rot;

	//a small adjustement: take the best rot from the IPs
	rottmp=fabs(angleCompute(p->uns[0],p->stb[0]));
	//printf("test1: %i %f\n",p->rot,rottmp);
	if (fabs(angleCompute(p->next->uns[0],p->next->stb[0]))<rottmp) {
		//printf("test2: %i %f\n",p->rot,rottmp);
		//printf("PIP switch requested on this Lobe...\n");
		cur=l;
		while (cur!=NULL) {
			cur->rot=-p->next->rot;
			cur=cur->next;
		}
	}


	return(l);
}

struct pip *getPips(struct chain *u, struct chain *s, int *np) {
	struct pip *p=NULL;
	struct pip *currentpip=NULL;
	struct chain *cu, *cs;
	struct chain *ends=NULL;
	char lmode[2048];
	cu=u;
	*np=0;

	printf("getpips\n");

	ends=s;
	while (ends->next!=NULL) {
		ends=ends->next;
	}
	ends=ends->prev;

	while (cu->next!=NULL) {
		cs=ends;
		while (cs!=NULL) {
			if (chainInter(cu,cs)) {
				//create a new pip
				if (currentpip==NULL) {
					p=(struct pip *)malloc(sizeof(struct pip));
					currentpip=p;
				} else {
					currentpip->next=(struct pip *)malloc(sizeof(struct pip));
					currentpip=currentpip->next;
				}
				currentpip->next=NULL;
				//add the info
				*np=*np+1;
				currentpip->uns[0]=cu;
				currentpip->uns[1]=cu->next;
				currentpip->stb[0]=cs;
				currentpip->stb[1]=cs->next;
				interCompute(cu,cs,&(currentpip->x[0]),&(currentpip->x[1]));
				rotCompute(cu,cs,&(currentpip->rot));
				sprintf(lmode,"pip at %s %s",STDOUT_FORMAT,STDOUT_FORMAT);
				printf(lmode,currentpip->x[0],currentpip->x[1]);
				printf(" rot=%i\n",currentpip->rot);
//				rotCompute(cu,cs,&(currentpip->rot));
//				printf(" or %i\n",currentpip->rot);
				//restrict stable manifold
				//ends=cs->prev;
				//changed to this on Fri Mar 12 for test
				//ends=cs->next;
				ends=cs;
			}
			cs=cs->prev;
		}
		cu=cu->next;
	}
	printf("getpips end\n");
	return(p);
}

double angleCompute(struct chain *c1, struct chain *c2) {
	double v1[2], v2[2];
	double rd;
	double n1, n2;
	v1[0]=c1->next->x[0]-c1->x[0];
	v1[1]=c1->next->x[1]-c1->x[1];
	v2[0]=c2->next->x[0]-c2->x[0];
	v2[1]=c2->next->x[1]-c2->x[1];
	rd=v1[0]*v2[0]+v1[1]*v2[1];
	n1=v1[0]*v1[0]+v1[1]*v1[1];
	n2=v2[0]*v2[0]+v2[1]*v2[1];
	rd=(n1>((double)(0.0)))?rd/sqrt(n1):rd;
	rd=(n2>((double)(0.0)))?rd/sqrt(n2):rd;
	return(rd);
}

void rotCompute(struct chain *c1, struct chain *c2, int *r) {
	double v1[2], v2[2];
	double rd;
	v1[0]=c1->next->x[0]-c1->x[0];
	v1[1]=c1->next->x[1]-c1->x[1];
	v2[0]=c2->next->x[0]-c2->x[0];
	v2[1]=c2->next->x[1]-c2->x[1];
	rd=v1[0]*v2[1]-v1[1]*v2[0];
	if (rd<0.0) {
		*r=-1;
	} else {
		*r=1;
	}
}

void interCompute(struct chain *c1, struct chain *c2, double *x, double *y) {
	/***check is c1 intersect c2***/
	double a1,b1,cc1;
	double a2,b2,cc2;
	double dx, dy;
	double det;

	//first line
	dx=c1->x[0]-c1->next->x[0];
	dy=c1->x[1]-c1->next->x[1];

	if ((fabs(dx)<=0.0)&&(fabs(dy)<=0.0)) {
	  printf("hoho!\n");
	}

	//need eqn of line c1
	if (fabs(dx)>=fabs(dy)) {
	  b1=1;
	  a1=-dy/dx;
	} else {
	  a1=1;
	  b1=-dx/dy;
	}
	cc1=-a1*c1->x[0]-b1*c1->x[1];

	//second line
	dx=c2->x[0]-c2->next->x[0];
	dy=c2->x[1]-c2->next->x[1];

	if ((fabs(dx)<=0.0)&&(fabs(dy)<=0.0)) {
	  printf("hoho!\n");
	}

	if (fabs(dx)>=fabs(dy)) {
		b2=1;
		a2=-dy/dx;
	} else {
		a2=1;
		b2=-dx/dy;
	}
	cc2=-a2*c2->x[0]-b2*c2->x[1];

	//intersection
	det=a1*b2-a2*b1;
	//printf("%15.10lf\t%15.10lf\n",a1,b2);
	//printf("%15.10lf\t%15.10lf\n",a2,b1);
	//printf("%15.10lf\n",det);
	if (fabs(det)>0.0) {
	  *x=cc1*b2-cc2*b1;
	  *x/=-det;
	  *y=a1*cc2-a2*cc1;
	  *y/=-det;
	} else {
	  printf("det=%e...fixed!\n",det);
	  *x=c1->x[0]+c1->next->x[0];
	  *x=*x/2.0;
	  *y=c1->x[1]+c1->next->x[1];
	  *y=*y/2.0;
	}
}

int chainInter(struct chain *c1, struct chain *c2) {
	/***check is c1 intersect c2***/
	double a,b,c;
	double dx, dy;
	double p1, p2;
	dx=c1->x[0]-c1->next->x[0];
	dy=c1->x[1]-c1->next->x[1];

	if ((fabs(dx)<=0.0)&&(fabs(dy)<=0.0)) {
	  return(0);
	}

	/*
	if ((fabs(dx)<1e-4)&&(fabs(dy)<1e-4)) {
	  return(0);
	}
	*/

	//need eqn of line c1
	if (fabs(dx)>=fabs(dy)) {
		b=1;
		a=-dy/dx;
	} else {
		a=1;
		b=-dx/dy;
	}
	c=-a*c1->x[0]-b*c1->x[1];
	//product for c2
	p1=a*c2->x[0]+b*c2->x[1]+c;
	p2=a*c2->next->x[0]+b*c2->next->x[1]+c;
	//test 1
	if (p1*p2>=0.0) {
		return(0);
	}
	//need eqn of line c2
	dx=c2->x[0]-c2->next->x[0];
	dy=c2->x[1]-c2->next->x[1];

	if ((fabs(dx)<=0.0)&&(fabs(dy)<=0.0)) {
	  return(0);
	}

	if (fabs(dx)>=fabs(dy)) {
		b=1;
		a=-dy/dx;
	} else {
		a=1;
		b=-dx/dy;
	}
	c=-a*c2->x[0]-b*c2->x[1];
	//product for c1
	p1=a*c1->x[0]+b*c1->x[1]+c;
	p2=a*c1->next->x[0]+b*c1->next->x[1]+c;
	//test 1
	if (p1*p2>=0.0) {
		return(0);
	}
	return(1);
}

void outputMani(FILE *f, struct chain *p) {
	/***output the content of the chain to a file ***/
	/***theEnd = end of manifold ; all manifold =? theEnd=NULL***/
	struct chain *cur;
	long n=0;
	char omode[2048];
	cur=p;
	sprintf(omode,"%s\t%s\n",OUT_FORMAT,OUT_FORMAT);
	while (cur!=NULL) {
		n++;
		fprintf(f,omode,cur->x[0],cur->x[1]);
		cur=cur->next;
	}
	printf("%i points\n",n);
}

void outputLobe(FILE *f, struct lobe *p) {
	/***output the content of the chain to a file ***/
	struct lobe *cur;
	char omode[2048];
	cur=p;
	sprintf(omode,"%s\t%s\t1.0\n",OUT_FORMAT,OUT_FORMAT);
	while (cur!=NULL) {
		fprintf(f,omode,cur->x[0],cur->x[1]);
		cur=cur->next;
	}
	fprintf(f,"0.0\t0.0\t0.0\n");
}

void outputPips(FILE *f, struct pip *p) {
	/***output the content of the pip struct to a file ***/
	struct pip *cur;
	char omode[2048];
	cur=p;
	sprintf(omode,"%s\t%s\n",OUT_FORMAT,OUT_FORMAT);
	while (cur!=NULL) {
		fprintf(f,omode,cur->x[0],cur->x[1]);
		cur=cur->next;
	}
}

long loadManiData(FILE *fil, struct chain **p) {
	//load points in the structure
	long npts=0;
	char *line=NULL;
	struct chain *cur=NULL;
	struct chain *pre=NULL;
	int i,j,k;
	double f;
	int stop;
	fpos_t pr;
	int errSet;
	int ret;

	//printf("  ***reading zone header\n");
	if (getLine(fil,&line)) {
		printf("  ***error generated by getLine\n");
		exit(0);
	}
	printf("  ***ZoneHeader=%sEOL\n",line);
	if (TecZone2IJ(&nx,&ny,line)<0) {
		printf("  ***error generated by TecZone2IJ\n");
		exit(0);
	}
	printf("  ***Zone size=%ix%i --> ",nx,ny);
	npts=nx;
	nx=-1;
	errSet=0;
	for (j=0;j<ny;j++) {
	  i=0;
	  stop=0;
	  while (stop==0) {
	    //		for (i=0;i<nx;i++) {
	    //add a new point
	    if (cur==NULL) {
	      *p=(struct chain *)malloc(sizeof(struct chain));
	      cur=*p;
	      (*p)->prev=NULL;
	    } else {
	      cur->next=(struct chain *)malloc(sizeof(struct chain));
	      cur->next->prev=cur;
	      cur=cur->next;
	    }
	    cur->next=NULL;
	    fgetpos(fil,&pr);
	    for (k=0;k<nVar;k++) {
	      ret=fscanf(fil,"%lf",&f);
	      //printf("%100.60lf\t",f);
	      if (ret>0) {
		cur->x[k]=f;
	      } else {
		if (errSet==0) {
		  //need to delete the last point (why?)
		  pre=cur->prev;
		  free(cur);
		  pre->next=NULL;
		  cur=pre;
		  i--;
		}
		errSet=1;
		stop=1;
	      }
	    }
	    //printf("\n");
	    i=i+1;
	    //getLine(fil,&line);
	    if (nx>0) {
	      if (i>=nx) {
		stop=1;
	      }
	    } else {
	      if (feof(fil)!=0) {
		if (stop==0) {
		  //need to delete the last point (why?)
		  pre=cur->prev;
		  free(cur);
		  pre->next=NULL;
		  cur=pre;
		  i--;
		}
		stop=1;
	      }
	    }
	  }
	  npts=i;
	  nx=npts;
	}
	if (errSet!=0) {
	  fsetpos(fil,&pr);
	}
	getLine(fil,&line);
	printf("%ix%i\n",nx,ny);
	free(line);
	return(npts);
}

void eraseChain(struct chain *p) {
	struct chain *nxt;
	struct chain *cur;
	cur=p;
	while (cur!=NULL) {
		nxt=cur->next;
		free(cur);
		cur=nxt;
	}
}

void erasePips(struct pip *p) {
	struct pip *nxt;
	struct pip *cur;
	cur=p;
	while (cur!=NULL) {
		nxt=cur->next;
		free(cur);
		cur=nxt;
	}
}

void eraseLobe(struct lobe *p) {
	struct lobe *nxt;
	struct lobe *cur;
	cur=p;
	while (cur!=NULL) {
		nxt=cur->next;
		free(cur);
		cur=nxt;
	}
}

FILE *setunsmanifold(const char *fn) {
	/***put the unstable manifold at the correct loc***/
	long z;
	int i,j,k;
	FILE *fil;
	char *line=NULL;
	double f;
	fpos_t p;
	int errSet;
	int ret;

	fil=fopen(fn,"r");
	if (fil==NULL) {
		printf("cannot open %s\n",fn);
		exit(0);
	}

	printf("  ***reading title\n");
	if (getLine(fil,&line)) {
		printf("  ***error generated by getLine\n");
		exit(0);
	}
	
	printf("  ***DataSetTitle=%sEOL\n",line);
	nVar=TecTitle2nVar(line);
	//free(line);
	if (nVar<0) {
	printf("  ***error generated by TecTitle2nVar\n");
		exit(0);
	}
	printf("  ***there are %i variables loaded\n",nVar);

	for (z=0;z<unsskip;z++) {
		//read header
		//printf("  ***reading zone header\n");
		if (getLine(fil,&line)) {
			printf("  ***error generated by getLine\n");
			exit(0);
		}
		printf("  ***ZoneHeader=%sEOL\n",line);
		if (TecZone2IJ(&nx,&ny,line)<0) {
			printf("  ***error generated by TecZone2IJ\n");
			exit(0);
		}
		//free(line);
		printf("  ***Zone size=%ix%i --> ",nx,ny);
		//read zone
		/*
		for (j=0;j<ny;j++) {
			for (i=0;i<nx;i++) {
				for (k=0;k<nVar;k++) {
					fscanf(fil,"%lf",&f);
				}	
			}
		}
		getLine(fil,&line);
		*/

		errSet=0;
		for (j=0;j<ny;j++) {
		  i=0;
		  while ((i<nx)||(nx<0)) {
		    //for (i=0;i<nx;i++) {
		    fgetpos(fil,&p);
		    for (k=0;k<nVar;k++) {
		      ret=fscanf(fil,"%lf",&f);
		      if (ret<=0) {
			nx=i;
			errSet=1;
		      }
		    }
		    i++;
		  }
		}
		if (errSet!=0) {
		  fsetpos(fil,&p);
		}
		getLine(fil,&line);

		printf("%ix%i\n",nx,ny);

	}
	free(line);
	return(fil);
}

void markstbmanifold(long znes, long *mark, const char *fn) {
	/***MARK the position of the requested zones in the array***/
	long z;
	int i,j,k;
	FILE *fil;
	char *line=NULL;
	double f;
	long position;
	long mem=0;
	int ret;
	fpos_t p;
	int errSet=0;

	fil=fopen(fn,"r");
	if (fil==NULL) {
		printf("cannot open %s\n",fn);
		exit(0);
	}

	printf("  ***reading title\n");
	if (getLine(fil,&line)) {
		printf("  ***error generated by getLine\n");
		exit(0);
	}
	printf("  ***DataSetTitle=%sEOL\n",line);
	nVar=TecTitle2nVar(line);
	//free(line);
	if (nVar<0) {
	  printf("  ***error generated by TecTitle2nVar\n");
	  exit(0);
	}
	printf("  ***there are %i variables loaded\n",nVar);

	for (z=0;z<stbskip;z++) {
		//read header
		//printf("  ***reading zone header\n");
		if (getLine(fil,&line)) {
			printf("  ***error generated by getLine\n");
			exit(0);
		}
		printf("  ***ZoneHeader=%sEOL\n",line);
		if (TecZone2IJ(&nx,&ny,line)<0) {
			printf("  ***error generated by TecZone2IJ\n");
			exit(0);
		}
		//free(line);
		printf("  ***Zone size=%ix%i --> ",nx,ny);

		//read zone

		errSet=0;
		for (j=0;j<ny;j++) {
		  i=0;
		  while ((i<nx)||(nx<0)) {
		    //for (i=0;i<nx;i++) {
		    fgetpos(fil,&p);
		    for (k=0;k<nVar;k++) {
		      ret=fscanf(fil,"%lf",&f);
		      if (ret<=0) {
			nx=i;
			errSet=1;
		      }
		    }
		    i++;
		  }
		}
		if (errSet!=0) {
		  fsetpos(fil,&p);
		}
		getLine(fil,&line);

		printf("%ix%i\n",nx,ny);

	}
	printf("marking manifold position file...\n");
	for (z=0;z<znes;z++) {
		//mark position
		position=ftell(fil);
		if (position<mem) {
			printf("position = %i memory = %i !\n",position,mem);
			exit(0);
		}
		mem=position;
		//fgetpos(fil,&(mark[znes-z+1]));
		mark[znes-z-1]=position;
		//read header
		//printf("  ***reading zone header\n");
		if (getLine(fil,&line)) {
			printf("  ***error generated by getLine\n");
			exit(0);
		}
		printf("  ***ZoneHeader=%sEOL\n",line);
		if (TecZone2IJ(&nx,&ny,line)<0) {
			printf("  ***error generated by TecZone2IJ\n");
			exit(0);
		}
		//free(line);
		printf("  ***Zone size=%ix%i --> ",nx,ny);
		//mark position
		//printf("reading position\n");
		//printf("done\n");
		//read zone
		/*
		for (j=0;j<ny;j++) {
			for (i=0;i<nx;i++) {
				for (k=0;k<nVar;k++) {
					fscanf(fil,"%lf",&f);
				}	
			}
		}
		//printf("yo!\n");
		getLine(fil,&line);
		//printf("looping...\n");
		*/

		errSet=0;
		for (j=0;j<ny;j++) {
		  i=0;
		  while ((i<nx)||(nx<0)) {
		    //for (i=0;i<nx;i++) {
		    fgetpos(fil,&p);
		    for (k=0;k<nVar;k++) {
		      ret=fscanf(fil,"%lf",&f);
		      if (ret<=0) {
			nx=i;
			errSet=1;
		      }
		    }
		    i++;
		  }
		}
		if (errSet!=0) {
		  fsetpos(fil,&p);
		}
		getLine(fil,&line);

		printf("%ix%i\n",nx,ny);


	}

	free(line);
	fclose(fil);
}

#define CURS_LNG 60
int CURS_POS=0;

void fancy_display_reset(const char *txt) {
  fprintf(stderr,"%s|",txt);
  fflush(stderr);
  CURS_POS=0;
}

void fancy_display_advance(long i, long n) {
  int w;
  int q;
  w=(int)(CURS_LNG*(((double)(i))/((double)(n))));
  w=(w>CURS_LNG)?CURS_LNG:w;
  for (q=CURS_POS+1;q<=w;q++) {
    fprintf(stderr,"=");
    fflush(stderr);
  }
  CURS_POS=w;
}

void fancy_display_stop(void) {
  int q;
    for (q=CURS_POS+1;q<=CURS_LNG;q++) {
    fprintf(stderr,"=");
    fflush(stderr);
  }
  CURS_POS=CURS_LNG;
  fprintf(stderr,"|\n");
  fflush(stderr);
}

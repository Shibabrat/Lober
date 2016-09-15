
#include <math.h>
#include "inandout.h"
#include <isInside.h>
#include <stdlib.h>

double *ar=NULL;

signed_curve *chain2signedcurve(struct chain *lob) {
  /*** This converts the chain* structure to a reg signed_curve from
       libisInside ***/
  long npt=0;
  struct chain *cur;
  signed_curve *c;
  long i;

  cur=lob;
  while (cur->next!=NULL) {
    npt++;
    ar=(double *)realloc(ar,2*npt*sizeof(double));
    ar[2*npt-2]=cur->x[0];
    ar[2*npt-1]=cur->x[1];
    cur=cur->next;
  }
  
  c=new signed_curve(npt,ar);
  return(c);
}

signed_curve *lobe2signedcurve(struct lobe *lob) {
  /*** This converts the chain* structure to a reg signed_curve from
       libisInside ***/
  long npt=0;
  struct lobe *cur;
  signed_curve *c;
  long i;

  cur=lob;
  while (cur->next!=NULL) {
    npt++;
    ar=(double *)realloc(ar,2*npt*sizeof(double));
    ar[2*npt-2]=cur->x[0];
    ar[2*npt-1]=cur->x[1];
    cur=cur->next;
  }
  
  c=new signed_curve(npt,ar);
  return(c);
}

/**************************************************************************
int isInside_wrapper_chain(int rot, double x, double y, struct chain *lob) {
  long npt=0;
  struct chain *cur;
  //double **ar=NULL;
  signed_curve *c;
  int ret;
  double pt[2];
  long i;

  cur=lob;
  while (cur->next!=NULL) {
    npt++;
    ar=(double *)realloc(ar,2*npt*sizeof(double));
    ar[2*npt-2]=cur->x[0];
    ar[2*npt-1]=cur->x[1];
    cur=cur->next;
  }

  c=new signed_curve(npt,ar);
  pt[0]=x;
  pt[1]=y;
  ret=c->isInside(pt);
  delete c;
  for (i=0;i<npt;i++) {
    //delete(ar[i]);
  }
  //delete(ar);
  return(ret);
}


int isInside_wrapper_lobe(double x, double y, struct lobe *lob) {
  long npt=0;
  struct lobe *cur;
  //double **ar=NULL;
  signed_curve *c;
  int ret;
  double pt[2];
  long i;

  cur=lob;
  while (cur->next!=NULL) {
    npt++;
    ar=(double *)realloc(ar,2*npt*sizeof(double));
    ar[2*npt-2]=cur->x[0];
    ar[2*npt-1]=cur->x[1];
    cur=cur->next;
  }

  c=new signed_curve(npt,ar);
  pt[0]=x;
  pt[1]=y;
  ret=c->isInside(pt);
  delete c;
  for (i=0;i<npt;i++) {
    //delete(ar[i]);
  }
  //delete(ar);
  return(ret);
}
********************************************************************/

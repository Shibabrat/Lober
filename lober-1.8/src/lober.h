#ifndef LOBER_LOBER_H
#define LOBER_LOBER_H
#include <stdio.h>

extern int nVar;

struct chain {
	double x[2];
	struct chain *next;
	struct chain *prev;
};

struct pip {
	struct chain *stb[2];
	struct chain *uns[2];
	struct pip *prev;
	struct pip *next;
	double x[2];
	int rot;
};

struct lobe {
	double x[2];
	struct lobe *next;
	int rot;
};

struct lobe *getSepar(struct chain *u, struct chain *s, struct pip *bip);
void outputSepar(struct lobe *sep);
void getLobeStat(struct lobe *l);
void getLobeStat(struct lobe *l, struct pip *p1, struct pip *p2, struct lobe *sep, struct pip *bip);
void markstbmanifold(long znes, long *mark, const char *fn);
FILE *setunsmanifold(const char *fn);
void eraseChain(struct chain *p);
void erasePips(struct pip *p);
void eraseLobe(struct lobe *p);
long loadManiData(FILE *f, struct chain **p);
void outputMani(FILE *f, struct chain *p);
void outputPips(FILE *f, struct pip *p);
void outputLobe(FILE *f, struct lobe *p);
void outputLobeb(FILE *f, long zonen, int nlobes, struct lobe **p);
struct pip *getPips(struct chain *u, struct chain *s, int *np);
int chainInter(struct chain *c1, struct chain *c2);
void interCompute(struct chain *c1, struct chain *c2, double *x, double *y);
void rotCompute(struct chain *c1, struct chain *c2, int *r);
struct lobe **getLobes(int npips, struct pip *pips, int *nlobes);
struct lobe *getLobe(struct pip *p);
//int isInLobes(double cx, double cy, int nlobes, struct lobe **l);
//int isInLobe(double cx, double cy, struct lobe *l);
int isInLobes(double cx, double cy, int nlobes, class signed_curve **l);
int isInLobe(double cx, double cy, class signed_curve *l);
void decimateLobe(double dmin, struct lobe *l);
long decimateMani(double dmin, struct chain *l, struct chain *theEnd);
void addMani(double x, double y, struct chain **c);
double angleCompute(struct chain *c1, struct chain *c2);

void fancy_display_reset(const char *txt);
void fancy_display_advance(long i, long n);
void fancy_display_stop(void);

#endif

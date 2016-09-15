#include <stdio.h>
#define getLineSegmLng 10
#define getLineEOL '\n'
int getLine(FILE *file, char **lne);
int strNFind(int shift,char *s1, char *s2);
int TecTitle2nVar(char *title);
int TecZone2IJ(int *nx, int *ny, char *title);
int setUpData(FILE *in, int *nx, int *ny, int *nVar, double ***d);
int strFind(char *s1, char *s2);
int strNotNFind(int shift,char *s1, char s2);

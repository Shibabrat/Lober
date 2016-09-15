#include <stdlib.h>
#include <string.h>

#define dfltSTR "unknown"

char *dflt(void) {
  char *c;
  int n;
  n=strlen(dfltSTR);
  n+=3;
  c=(char *)malloc(n*sizeof(char));
  strcpy(c,dfltSTR);
  return(c);
}

char *getStringCName(void) {
  char *c;
  c=getenv("HOST");
  if (c==NULL) {
    c=getenv("HOSTNAME");
  }
  if (c==NULL) {
    c=dflt();
  }
  return(c);
}

char *getStringOSName(void) {
  char *c;
  char *mach;
  char *str;
  c=getenv("OS");
  if (c==NULL) {
    c=getenv("OSTYPE");
  }

  if (c!=NULL) {
    mach=getenv("MACH");
    if (mach==NULL) {
      mach=getenv("MACHTYPE");
    }
    if (mach!=NULL) {
      str=(char *)malloc((strlen(c)+strlen(mach)+2)*sizeof(char));
      strcpy(str,c);
      strcat(str,"-");
      strcat(str,mach);
      c=str;
    }
  }

  if (c==NULL) {
    c=dflt();
  }

  return(c);
}

char *getUserName(void) {
  char *c;
  c=getenv("USER");
  if (c==NULL) {
    c=getenv("USERNAME");
  }
  if (c==NULL) {
    c=dflt();
  }
  return(c);
}

char *getUserGroup(void) {
  char *c;
  c=getenv("GROUP");
  if (c==NULL) {
    c=getenv("GROUPNAME");
  }
  if (c==NULL) {
    c=dflt();
  }
  return(c);
}

char *getPath(void) {
  char *c;
  c=getenv("PWD");
  if (c==NULL) {
    c=dflt();
  }
  return(c);
}

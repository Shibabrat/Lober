#include <isInside.h>
#include "version.h"
char *       __IS_INSIDE__NAME__(void) {
  return(__LIBNAME__);
}

unsigned int __IS_INSIDE__VERSION_INTERFACE__(void) {
  return(__LIBINTERFACE__);
}

unsigned int __IS_INSIDE__VERSION_REVISION__(void) {
  return(__LIBREVISION__);
}

unsigned int __IS_INSIDE__VERSION_COMPATIBILITY__(void) {
  return(__LIBCOMPATIBILITY__);
}

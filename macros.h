#ifndef MACROS
#define MACROS

#include <stdio.h>
#include <stdlib.h>


#define NEW(objectptrvar) do { void *_tmp=malloc(sizeof(*(objectptrvar))); if(_tmp==NULL) { printf("Main memory allocation error of %lu bytes in %s:%d\n",(long unsigned)sizeof(*(objectptrvar)),__FILE__,__LINE__); abort(); } (objectptrvar)=_tmp; } while(0)

#define DELETE(objectptrvar) free(objectptrvar)

#define CPU 0
#define MPI 1


#endif

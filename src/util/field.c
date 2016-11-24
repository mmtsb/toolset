#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>

#include "field.h"


// ##### Field ######################################################

const int Field::maxfields=200;

int Field::contains(char c, char *set) {
  for (char *iptr=set; *iptr!=0; iptr++) {
    if (c==*iptr) {
      return 1;
    }
  }
  return 0;
}


Field::Field(char *inp, char *sep) : nfields(0), first(0) {
  tstr=new char[strlen(inp)+2];

  strcpy(tstr,inp);

  char *cptr=tstr;
  do {
    for (; *cptr!=0 && contains(*cptr,sep); cptr++);
    if (*cptr!=0) {
      fstr[nfields++]=cptr;
      for (;*cptr!=0 && !contains(*cptr,sep); cptr++);
      if (*cptr!=0) {
	*cptr=0;
	cptr++;
      }
//      fprintf(stderr,"nfields=%d, str=%s\n",nfields-1,fstr[nfields-1]);
    }
  } while (*cptr!=0);
}


Field::~Field() {
  delete tstr;
}

void Field::shift() {
  first++;
}

char *Field::operator[](int inx) {
  return fstr[inx-first];
}











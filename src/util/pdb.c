#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "vector.h"
#include "pdb.h"

// *** PDBEntry ******

PDBEntry::PDBEntry() : 
 atminx(0), resnum(0), alt(0), chn(0), coor(0.0,0.0,0.0), 
 aux1(0.0), aux2(0.0), selFlag(0) {
  *atmname=0;
  *resname=0;
  *segid=0;
}

PDBEntry::PDBEntry(int ainx, char *aname, char *rname, int rnum, 
		   Vector v, char *seg) :
 atminx(ainx), resnum(rnum), alt(0), chn(0), coor(v), 
 aux1(0.0), aux2(0.0), selFlag(0) {
  strcpy(atmname,aname);
  strcpy(resname,rname);
  if (seg!=0) {
    strcpy(segid,seg); 
  } else {
    *segid=0;
  }
}

PDBEntry& PDBEntry::operator=(const PDBEntry& from) {
  atminx=from.atminx;
  resnum=from.resnum;
  alt=from.alt;
  chn=from.chn;
  coor=from.coor;
  aux1=from.aux1;
  aux2=from.aux2;
  selFlag=from.selFlag;
  if (from.atmname!=0) strcpy(atmname,from.atmname);
  if (from.resname!=0) strcpy(resname,from.resname);
  if (from.segid!=0) strcpy(segid,from.segid);
  return (*this);
}

char PDBEntry::type() {
  char *cptr;
  for (cptr=atmname; *cptr!=0 && isdigit(*cptr); cptr++);
  return *cptr;
}

char *PDBEntry::substr(char *str, int offset, int len) {
  static char buffer[1024];
  
  char *cptr,*dptr;
  
  cptr=&str[offset];
  dptr=buffer;
  for (int i=0; i<len; i++) {
    if (*cptr!=' ') 
      *dptr++=*cptr;
    cptr++;
  }
  *dptr=0;

  return buffer;
}

int PDBEntry::read(FILE *fptr, SelEnum selMode) {
  char line[1024];

  if (!fgets(line,1024,fptr))
    return -1;

  if (strlen(line)<50 || (strncmp(line,"ATOM",4) && strncmp(line,"HETATM",6)))
    return 0;

  atminx=atoi(substr(line,6,6));
  alt=line[16];
  if (isdigit(alt)) {
    strcpy(atmname,substr(line,12,5));
  } else {
    strcpy(atmname,substr(line,12,4));
  }
  strcpy(resname,substr(line,17,4));
  resnum=atoi(substr(line,22,5));
  chn=line[21];

  coor=Vector(atof(substr(line,30,8)),
  	      atof(substr(line,38,8)),
  	      atof(substr(line,46,8)));
  aux1=atof(substr(line,54,6));
  aux2=atof(substr(line,60,6));
  strcpy(segid,substr(line,72,4));

  if (selMode==ALL ||
      (type()!='H' && (selMode==ALL || selMode==HEAVY)) ||
      (!strcmp(atmname,"CA") && (selMode!=CB || !strcmp(resname,"GLY"))) ||
      (!strcmp(atmname,"CB") && selMode!=CA)) 
    selFlag=1;
  else 
    selFlag=0;

  return 1;
}

void PDBEntry::write(FILE *fptr) {
  if (strlen(atmname)>3) {
    if (resnum>99999) {
     fprintf(fptr,
      "ATOM%7d %-4s%c%-4s%c%6d  %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",
      atminx,atmname,(alt==0)?' ':alt,resname,(chn==0)?' ':chn,resnum,coor.x(),coor.y(),coor.z(),
      aux1,aux2,segid);
    } else {
     fprintf(fptr,
      "ATOM%7d %-4s%c%-4s%c%5d   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",
      atminx,atmname,(alt==0)?' ':alt,resname,(chn==0)?' ':chn,resnum,coor.x(),coor.y(),coor.z(),
      aux1,aux2,segid);
    }
  } else {
    if (resnum>99999) {
     fprintf(fptr,
      "ATOM%7d  %-3s%c%-4s%c%6d  %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",
      atminx,atmname,(alt==0)?' ':alt,resname,(chn==0)?' ':chn,resnum,coor.x(),coor.y(),coor.z(),
      aux1,aux2,segid);
    } else {
     fprintf(fptr,
      "ATOM%7d  %-3s%c%-4s%c%5d   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",
      atminx,atmname,(alt==0)?' ':alt,resname,(chn==0)?' ':chn,resnum,coor.x(),coor.y(),coor.z(),
      aux1,aux2,segid);
    }

  }
}


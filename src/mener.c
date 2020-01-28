// 2002, Michael Feig, Brooks group, TSRI
// SICHO model lattice energy

// *** C++ code ***
//
// requires:
//   field.h field.c
//   melib.f
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "field.h"

extern "C" {
#if defined IBM
  void msetup(char *ddir, int *lenf, int *seq, int *sec,
	       double *aarep, double *compr, double *esc0,
	       double *arl0, double *enone, double *ehbn,
	       double *ars, double *burek, double *enscal,
	       double *es3, double *temp, int ddirlen);
  void mener(int *cx, int *cy, int *cz, double *ener, double *eshort,
	      double *epair, double *eburial, double *ecorrect, double *ecent);
#else
  void msetup_(char *ddir, int *lenf, int *seq, int *sec,
	       double *aarep, double *compr, double *esc0,
	       double *arl0, double *enone, double *ehbn,
	       double *ars, double *burek, double *enscal,
	       double *es3, double *temp, int ddirlen);
  void mener_(int *cx, int *cy, int *cz, double *ener, double *eshort,
	      double *epair, double *eburial, double *ecorrect, double *ecent);
#endif
};

class SeqNum {
private:         
  static char sng[];
  static char nm[][4];
  static char mnm[][4];

public:
  static int fromSingle(char c);
  static int fromName(char *s);
  static int fromMName(char *s);

  static const char *name(int inx);
  static char single(int inx); 
};

char SeqNum::sng[]="ARNDCQEGHILKMFPSTWYV";
char SeqNum::nm[][4]={ "ALA", "ARG", "ASN", "ASP", "CYS",
		       "GLN", "GLU", "GLY", "HIS", "ILE",
		       "LEU", "LYS", "MET", "PHE", "PRO",
		       "SER", "THR", "TRP", "TYR", "VAL" };

char SeqNum::mnm[][4]={ "GLY", "ALA", "SER", "CYS", "VAL",
		        "THR", "ILE", "PRO", "MET", "ASP",
		        "ASN", "LEU", "LYS", "GLU", "GLN",
		        "ARG", "HIS", "PHE", "TYR", "TRP" };

int SeqNum::fromSingle(char c) {
  for (int i=0; i<20; i++) {
    if (sng[i] == c)
      return i;
  }
  return -1;
}

int SeqNum::fromName(char *s) {
  for (int i=0; i<20; i++) {
    if (!strcmp(nm[i],s)) 
      return i;
  }
  return -1;
}

int SeqNum::fromMName(char *s) {
  for (int i=0; i<20; i++) {
    if (!strcmp(mnm[i],s))
      return i;
  }
  return -1;
}
 
char SeqNum::single(int inx) {
  return (inx>=0 && inx<20)?sng[inx]:'\0';
}

const char *SeqNum::name(int inx) {
  return (inx>=0 && inx<20)?nm[inx]:"";
}

class Sequence {
private:
  int *res;
  int *mres;
  int *sec;
  int n;
  int maxn;

public:
  Sequence(int max);
  ~Sequence();

  void readMONSSTERFile();

  int residue(int inx);
  int secondary(int inx);
  int nres() { return n; }
  int *mresAddr() { return mres; }
  int *secAddr() { return sec;}
};

Sequence::Sequence(int max) : n(0), maxn(max) {
  res=new int[maxn];
  mres=new int[maxn];
  sec=new int[maxn];
}

Sequence::~Sequence() {
  delete res;
  delete mres;
  delete sec;
}

void Sequence::readMONSSTERFile() {
  char line[512];
  fgets(line,512,stdin);
  while (n<maxn) {
    fgets(line,512,stdin);
    Field f(line);
    res[n]=SeqNum::fromName(f[1]);
    mres[n]=SeqNum::fromMName(f[1]);
    sec[n]=atoi(f[2]);
    n++;
  }
}

int Sequence::residue(int inx) {
  return (inx>=0 && inx<n)?res[inx]:0;
}

int Sequence::secondary(int inx) {
  return (inx>=0 && inx<n)?sec[inx]:0;
}
class Chain {
 public:
  int len;
  int *x,*y,*z;

  Chain(int s) : len(s) {
    x=new int[s];
    y=new int[s];
    z=new int[s];
  }
  ~Chain() {
    delete x;
    delete y;
    delete z;
  }

  int idiff(int inx, int jnx);
  int colinear(int inx, int jnx, int knx);

  Chain& operator=(const Chain& from);

  void save(char *fname);
};

void Chain::save(char *fname) {
  FILE *fptr=fopen(fname,"w");
  fprintf(fptr,"%6d\n",len);
  for (int i=0; i<len; i++) {
    fprintf(fptr,"%6d%6d%6d\n",x[i],y[i],z[i]);
  }
  fclose(fptr);
}

Chain& Chain::operator=(const Chain& from) {
  for (int i=0; i<len; i++) {
    x[i]=from.x[i];
    y[i]=from.y[i];
    z[i]=from.z[i];
  }
  return (*this);
}

int Chain::idiff(int inx, int jnx) {
  int dx=x[inx]-x[jnx];
  int dy=y[inx]-y[jnx];
  int dz=z[inx]-z[jnx];
  return (dx*dx+dy*dy+dz*dz);
}

int Chain::colinear(int inx, int jnx, int knx) {
  int ix=x[inx]-x[jnx];
  int iy=y[inx]-y[jnx];
  int iz=z[inx]-z[jnx];

  int jx=x[jnx]-x[knx];
  int jy=y[jnx]-y[knx];
  int jz=z[jnx]-z[knx];

  int kx=iy*jz-iz*jy;
  int ky=jx*iz-ix*jz;
  int kz=ix*jy-iy*jx;

  return ((kx*kx+ky*ky+kz*kz)==0)?1:0;
}

int main(int argc, char **argv) {
  char datdir[256];
  fgets(datdir,256,stdin);

  double softcore,central,stiff,pair;
  double kdcore,hbond,pshort,burial;
  double multibody,threebody,temp;

  scanf("%lf%lf%lf",&softcore,&central,&stiff);
  scanf("%lf%lf%lf",&pair,&kdcore,&hbond);
  scanf("%lf%lf%lf",&pshort,&burial,&multibody);
  scanf("%lf%lf",&threebody,&temp);

  int lenf;

  scanf("%d",&lenf);

  Sequence s(lenf);
  s.readMONSSTERFile();
  
#if defined IBM
  msetup(datdir,&lenf,s.mresAddr(),s.secAddr(),&softcore,&central,
          &stiff,&pair,&kdcore,&hbond,&pshort,&burial,&multibody,&threebody,&temp,
	  strlen(datdir));
#else
  msetup_(datdir,&lenf,s.mresAddr(),s.secAddr(),&softcore,&central,
          &stiff,&pair,&kdcore,&hbond,&pshort,&burial,&multibody,&threebody,&temp,
	  strlen(datdir));
#endif

  char line[256];

  fgets(line,256,stdin);
  int clen=atoi(line);

  if (clen==lenf+2) 
    fgets(line,256,stdin);

  Chain chain(lenf);
  for (int i=0; i<lenf; i++) {
    scanf("%d%d%d",&chain.x[i],&chain.y[i],&chain.z[i]);
  }

  double ener,ecent,eshort,epair,eburial,ecorrect;
#if defined IBM
  mener(chain.x,chain.y,chain.z,&ener,&eshort,&epair,&eburial,&ecorrect,&ecent);
#else
  mener_(chain.x,chain.y,chain.z,&ener,&eshort,&epair,&eburial,&ecorrect,&ecent);
#endif
  
  printf("%15.8lf %15.8lf %15.8lf %15.8lf %15.8lf %15.8lf\n",
	 ener,eshort,epair,eburial,ecorrect,ecent);
}
  

// http://mmtsb.scripps.edu/doc/gdt.html
// 2001, Michael Feig (meikel@scripps.edu)
// The Scripps Research Institute

// *** C++ code ***
//
// requires:
//   vector.h, pdb.c, pdb.h, field.h, field.c, charmmlsq.f
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "vector.h"
#include "field.h"
#include "pdb.h"

void usage() {
  fprintf(stderr,"usage:   gdt [-cutoff value] [-ts] refPDB [compPDB]\n");
  exit(1);
}

extern "C" {
#if defined(IRIX) || defined(LINUX) || defined(SUNOS) || defined(ALPHA)
  void frotu_(double *r, double *u);
#else
#if defined(CRAY)
  void FROTU(double *r, double *u);
#else
  void frotu(double *r, double *u);
#endif
#endif
};


int readPDB(char *pdbName, PDBEntry* cskel) {
  FILE *fptr;
  fptr= (!strcmp(pdbName,"-") || *pdbName==0) ? stdin : fopen(pdbName,"r");

  if (fptr==0) {
    fprintf(stderr,"cannot open PDB file %s\n",pdbName);
    exit(1);
  }

  int natom=0;

  while (!feof(fptr)) {
    if (cskel[natom].read(fptr,CA)>0) {
      if (cskel[natom].selected()) natom++;
    }
  }

  fclose(fptr);
  
  //    fprintf(stderr,"finished reading %d CA atoms from %s\n",natom,pdbName);

  return natom;
}

int matchResidue(PDBEntry *e1, int i1, PDBEntry *e2, int i2) {
  char *r1=e1[i1].residueName();
  char *r2=e2[i2].residueName();

  if (e1[i1].residueNumber()==e2[i2].residueNumber()) {
    if (!strcmp(r1,r2) || 
	((!strcmp(r1,"HIS") || !strcmp(r1,"HSD") || !strcmp(r1,"HSE")) &&
	 (!strcmp(r2,"HIS") || !strcmp(r2,"HSD") || !strcmp(r2,"HSE")))) {
      return 1;
    } else {
      fprintf(stderr,"residue %d name mismatch %s/%s\n",
	      e1[i1].residueNumber(),r1,r2);
      return 1;
    }
  } 
  return 0;
}

int fit(PDBEntry *cmp, PDBEntry *ref, int ntot, int *map, 
	int *list, int& nlist, double cutoff) {
  Vector refskel[1000];
  Vector cmpskel[1000];

  int i;
  for (i=0; i<nlist; i++) {
    cmpskel[i]=cmp[list[i]].coordinates();
    refskel[i]=ref[map[list[i]]].coordinates();
    //    fprintf(stderr,"%d ",list[i]);
  }
  //  fprintf(stderr,"\n");

  Vector tc1(0,0,0);
  Vector tc2(0,0,0);

  for (i=0; i<nlist; i++) {
    tc1+=refskel[i];
    tc2+=cmpskel[i];
    //    fprintf(stderr,"%f %f %f\n",refskel[i].x(),refskel[i].y(),refskel[i].z());
    //    fprintf(stderr,"%f %f %f\n",cmpskel[i].x(),cmpskel[i].y(),cmpskel[i].z());
  }
  tc1/=(double)nlist;
  tc2/=(double)nlist;

  double r[9];
  double u[9];
  for (i=0; i<9; i++) 
    r[i]=0.0;

  for (i=0; i<nlist; i++) {
    refskel[i]-=tc1;
    cmpskel[i]-=tc2;

    r[0]+=cmpskel[i].x()*refskel[i].x();
    r[1]+=cmpskel[i].x()*refskel[i].y();
    r[2]+=cmpskel[i].x()*refskel[i].z();

    r[3]+=cmpskel[i].y()*refskel[i].x();
    r[4]+=cmpskel[i].y()*refskel[i].y();
    r[5]+=cmpskel[i].y()*refskel[i].z();

    r[6]+=cmpskel[i].z()*refskel[i].x();
    r[7]+=cmpskel[i].z()*refskel[i].y();
    r[8]+=cmpskel[i].z()*refskel[i].z();
  }

#if defined(IRIX) || defined(LINUX) || defined(SUNOS) || defined(ALPHA)
  frotu_(r,u);
#else
#if defined(CRAY)
  FROTU(r,u);
#else
  frotu(r,u);
#endif
#endif
  
  int different=0;

  int plist=nlist;

  nlist=0;
  for (i=0; i<ntot; i++) {
    Vector v(cmp[i].coordinates());
    v-=tc2;
    Vector nv(u[0]*v.x()+u[3]*v.y()+u[6]*v.z(),
	      u[1]*v.x()+u[4]*v.y()+u[7]*v.z(),
	      u[2]*v.x()+u[5]*v.y()+u[8]*v.z());
    nv+=tc1;

    Vector dv=nv-ref[map[i]].coordinates();
    double dist=dv.norm();

    //       fprintf(stderr,"checking %d %lf (%lf)\n",i,dist,cutoff);
    
    if (dist<=cutoff) {
      if (list[nlist]!=i || nlist>=plist) different=1;
      list[nlist++]=i;
    }
  }
  for (i=0; i<nlist; i++) {
    //    fprintf(stderr,"%d ",list[i]);
  }
  //  fprintf(stderr,"\ndifferent: %d\n",different);

  return different;
}

int doCutoff(PDBEntry *cmpent, PDBEntry *refent, int crefent, int *map, 
	     double cutoff) {
  
  int maxnlist=0;

  int list[1000];
  int nlist=0;
  for (int d=1; d<crefent/3; d++) {
    for (int i=0; i<crefent-d-d; i++) {
      //      fprintf(stderr,"i: %d  d: %d\n",i,d);
      nlist=0;
      list[nlist++]=i;
      list[nlist++]=i+d;
      list[nlist++]=i+d+d;

      int different;
      int cnt=0;
      do {
	different=fit(cmpent,refent,crefent,map,list,nlist,cutoff);
	if (nlist>maxnlist)
	  maxnlist=nlist;
      } while (different && nlist>=3 && ++cnt<10);
    }
  }
  return maxnlist;
}



int main(int argc, char **argv) {
  if (argc<1) 
    usage();

  int i;

  char reffile[2048];
  *reffile=0;
  char cmpfile[2048];
  *cmpfile=0;

  double cutoff=5.0;
  int ts=0;

  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i],"-help")) {
      usage();
    } else if (!strcmp(argv[i],"-cutoff")) {
      cutoff=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-ts")) {
      ts=1;
    } else {
      strcpy(reffile,argv[i++]);
      if (i<argc)
	strcpy(cmpfile,argv[i]);
    }
  }

  PDBEntry refent[5000];
  PDBEntry cmpent[5000];
  
  int nrefent=readPDB(reffile,refent);
  int crefent=readPDB(cmpfile,cmpent);

  int *map=new int[crefent];

  int ri=0;
  for (int ci=0; ci<crefent; ci++) {
    while (!matchResidue(cmpent,ci,refent,ri) && ri<nrefent) ri++;
    if (ri>=nrefent) {
      fprintf(stderr,"cannot find match for %s%d:CA, no GDT will be calculated\n",
	      cmpent[ci].residueName(),cmpent[ci].residueNumber());
      exit(1);
    }
    map[ci]=ri;
  }    

  if (!ts) {
    int num=doCutoff(cmpent,refent,crefent,map,cutoff);
    printf("%d %lf\n",num,(double)num/(double)crefent*100.0);    
  } else {
    int num1=doCutoff(cmpent,refent,crefent,map,1.0);
    int num2=doCutoff(cmpent,refent,crefent,map,2.0);
    int num4=doCutoff(cmpent,refent,crefent,map,4.0);
    int num8=doCutoff(cmpent,refent,crefent,map,8.0);
    double gdtts=((double)num1+(double)num2+(double)num4+(double)num8)*25.0/
      (double)crefent;
    printf("%lf\n",gdtts);
  }
}




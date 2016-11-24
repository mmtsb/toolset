// 2007, Michael Feig, Michigan State University

// *** C++ code ***
//
// requires:
//   field.h, field.c,
//   aa.h aa.c
//   pdb.h, pdb.c
//   vector.h
//
//   brent.c, f1dim.c, frprmn.c linmin.c mnbrak.c nrutil.h nrutil.c  
//   (modified from Numerical Recipes Software)
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "field.h"
#include "vector.h"
#include "pdb.h"
#include "aa.h"

void usage() {
  fprintf(stderr,"usage:    rebprimo [options] [pdbfile]\n");
  fprintf(stderr,"options:  [-d datadir]\n");
  exit(1);
}

const int MAXATOM=1000000;

void readPDB(char *pdbName, PDBEntry *pdb, int& natom, int *resStart, int& nres) {
  FILE *fptr;

  if (!strcmp(pdbName,"-")) {
    fptr=stdin;
  } else {
    fptr=fopen(pdbName,"r");
    if (fptr==0) {
      fprintf(stderr,"Cannot open PDB file %s\n",pdbName);
      exit(1);
    }
  }

  int lastnum=-1;
  nres=0;
  natom=0;
  
  while (!feof(fptr)) {
    if (pdb[natom].read(fptr)>0) {
      if (pdb[natom].residueNumber()!=lastnum) {
	resStart[nres]=natom;
	nres++;
	lastnum=pdb[natom].residueNumber();
      }
      natom++;
    }
  }

  fclose(fptr);
  resStart[nres]=natom;
}

PDBEntry *findPDB(char* name, PDBEntry *mol, int start, int max) {
  int resnum;

  resnum=mol[start].residueNumber();
  for (int i=start; i<max && mol[i].residueNumber()==resnum; i++) {
    if (!strcmp(mol[i].atomName(),name)) {
      return &mol[i];
    }
  }
  fprintf(stderr,"cannot find atom %s (residue %s:%d) in PDB\n",name,mol[start].residueName(),resnum);
  exit(1);
  return 0;
}

double dihedral(Vector& v1, Vector& v2, Vector& v3, Vector& v4) {
  Vector d12=v1-v2;
  Vector d23=v2-v3;
  Vector d43=v4-v3;

  Vector p1=d12.cross(d23);
  p1/=p1.norm();

  Vector p2=d43.cross(d23);
  p2/=p2.norm();

  double dp12=p1*p2;
  double angle=acos(dp12);

  Vector p3=p1.cross(p2);
  double dp233=p3*d23;

  if (dp233>0.0) {
    angle=-angle;
  }
  return angle*180.0/M_PI;
}

int main(int argc, char **argv) {
  char datadir[256];
  char *mptr;
  if ((mptr=getenv("MMTSBDIR"))!=(char *)0) {
    strcpy(datadir,mptr);
    strcat(datadir,"/data");
  } else 
    strcpy(datadir,"./");

  char inputfile[256];
  strcpy(inputfile,"-");

  int i,j;

  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i],"-d")) {
      strcpy(datadir,argv[++i]);
    } else if (!strcmp(argv[i],"-help")) {
      usage();
    } else {
      strcpy(inputfile,argv[i]);
      break;
    }
  }

  PDBEntry *mol;
  mol=new PDBEntry[MAXATOM];
  
  int nres,natom;
  int *rstart;
  rstart=new int[MAXATOM/2];

  readPDB(inputfile,mol,natom,rstart,nres);

  fprintf(stderr,"read %d atoms, %d residues from %s\n",
	  natom,nres,inputfile);

  AA **res;
  res=new AA*[nres];

  // initial setup of atoms from PRIMO 
  for (i=0; i<nres; i++) {
    PDBEntry *m=&mol[rstart[i]];
    res[i]=getAA(m->residueName());
    res[i]->pdbnum=m->residueNumber();
    res[i]->ca=findPDB("CA",mol,rstart[i],natom)->coordinates();
    res[i]->n=findPDB("N",mol,rstart[i],natom)->coordinates();
    res[i]->c=findPDB("CO",mol,rstart[i],natom)->coordinates();

    if (!strcmp(res[i]->name(),"ALA")) 
      res[i]->set("CB",findPDB("SC1",mol,rstart[i],natom)->coordinates());

    if (!strcmp(res[i]->name(),"VAL")) 
      res[i]->set("CG1",findPDB("SC1",mol,rstart[i],natom)->coordinates());

    if (!strcmp(res[i]->name(),"LEU")) {
      res[i]->set("CG",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("CD1",findPDB("SC2",mol,rstart[i],natom)->coordinates());
    }
	
    if (!strcmp(res[i]->name(),"ILE")) {
      res[i]->set("CG1",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("CD",findPDB("SC2",mol,rstart[i],natom)->coordinates());
    }
	
    if (!strcmp(res[i]->name(),"PRO")) 
      res[i]->set("CG",findPDB("SC1",mol,rstart[i],natom)->coordinates());
	
    if (!strcmp(res[i]->name(),"CYS")) 
      res[i]->set("SG",findPDB("SC1",mol,rstart[i],natom)->coordinates());
	
    if (!strcmp(res[i]->name(),"MET")) {
      res[i]->set("CG",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("SD",findPDB("SC2",mol,rstart[i],natom)->coordinates());
      res[i]->set("CE",findPDB("SC3",mol,rstart[i],natom)->coordinates());
    }

    if (!strcmp(res[i]->name(),"PHE")) {
      res[i]->set("CG",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("CE1",findPDB("SC2",mol,rstart[i],natom)->coordinates());
      res[i]->set("CE2",findPDB("SC3",mol,rstart[i],natom)->coordinates());
    }

    if (!strcmp(res[i]->name(),"TYR")) {
      res[i]->set("CG",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("CE1",findPDB("SC2",mol,rstart[i],natom)->coordinates());
      res[i]->set("CE2",findPDB("SC3",mol,rstart[i],natom)->coordinates());
      res[i]->set("OH",findPDB("SC4",mol,rstart[i],natom)->coordinates());
    }

    if (!strcmp(res[i]->name(),"TRP")) {
      res[i]->set("CG",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("NE1",findPDB("SC2",mol,rstart[i],natom)->coordinates());
      res[i]->set("CE3",findPDB("SC3",mol,rstart[i],natom)->coordinates());
      res[i]->set("CZ2",findPDB("SC4",mol,rstart[i],natom)->coordinates());
    }

    if (!strcmp(res[i]->name(),"SER")) 
      res[i]->set("OG",findPDB("SC1",mol,rstart[i],natom)->coordinates());

    if (!strcmp(res[i]->name(),"THR")) {
      res[i]->set("OG1",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("CG2",findPDB("SC2",mol,rstart[i],natom)->coordinates());
    }

    if (!strcmp(res[i]->name(),"ASN")) {
      res[i]->set("CB",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("OD1",findPDB("SC2",mol,rstart[i],natom)->coordinates());
      res[i]->set("ND2",findPDB("SC3",mol,rstart[i],natom)->coordinates());
    }

    if (!strcmp(res[i]->name(),"GLN")) {
      res[i]->set("CG",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("OE1",findPDB("SC2",mol,rstart[i],natom)->coordinates());
      res[i]->set("NE2",findPDB("SC3",mol,rstart[i],natom)->coordinates());
    }

    if (!strcmp(res[i]->name(),"LYS")) {
      res[i]->set("CG",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("CD",findPDB("SC2",mol,rstart[i],natom)->coordinates());
      res[i]->set("CE",findPDB("SC3",mol,rstart[i],natom)->coordinates());
      res[i]->set("NZ",findPDB("SC4",mol,rstart[i],natom)->coordinates());
    }

    if (!strcmp(res[i]->name(),"ARG")) {
      res[i]->set("CG",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("CD",findPDB("SC2",mol,rstart[i],natom)->coordinates());
      res[i]->set("NE",findPDB("SC3",mol,rstart[i],natom)->coordinates());
      res[i]->set("NH1",findPDB("SC4",mol,rstart[i],natom)->coordinates());
      res[i]->set("NH2",findPDB("SC5",mol,rstart[i],natom)->coordinates());
    }

    if (!strcmp(res[i]->name(),"HIS") || !strcmp(res[i]->name(),"HSD") || 
	!strcmp(res[i]->name(),"HSE") || !strcmp(res[i]->name(),"HSP")) {
      res[i]->set("CG",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("ND1",findPDB("SC2",mol,rstart[i],natom)->coordinates());
      res[i]->set("NE2",findPDB("SC3",mol,rstart[i],natom)->coordinates());
    }

    if (!strcmp(res[i]->name(),"ASP")) {
      res[i]->set("CB",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("OD1",findPDB("SC2",mol,rstart[i],natom)->coordinates());
      res[i]->set("OD2",findPDB("SC3",mol,rstart[i],natom)->coordinates());
    }

    if (!strcmp(res[i]->name(),"GLU")) {
      res[i]->set("CG",findPDB("SC1",mol,rstart[i],natom)->coordinates());
      res[i]->set("OE1",findPDB("SC2",mol,rstart[i],natom)->coordinates());
      res[i]->set("OE2",findPDB("SC3",mol,rstart[i],natom)->coordinates());
    }
  }

  // fix backbone carbonyl atoms C and O
  double dcac=1.49;
  double dcn=1.345;
  double dco=1.23;
  double dcac2=dcac*dcac;
  double dcn2=dcn*dcn;
  for (i=0; i<nres; i++) {
    if (i<nres-1 && res[i]->pdbnum == res[i+1]->pdbnum-1) {
      Vector rcaco=res[i]->c-res[i]->ca;
      Vector rcan=res[i+1]->n-res[i]->ca;
      double dcaco=rcaco.norm();
      double dcan=rcan.norm();
      Vector ncaco=rcaco/dcaco;
      Vector ncan=rcan/dcan;
      double dcaco2=dcaco*dcaco;
      double dcan2=dcan*dcan;
      double a=(dcan2+dcac2-dcn2)/(2.0*dcan);
      double h2=dcac2-a*a;
      double h=sqrt(h2);
      double cosa=ncaco*ncan; 
      double sina=sqrt(1.0-cosa*cosa);
      double x=h/sina;   
      double ab=sqrt(x*x-h2);
      double aa=a-ab;
      Vector c=res[i]->ca+ncan*aa+ncaco*x;
      Vector rcco=res[i]->c-c;
      rcco/=rcco.norm();
      res[i]->c=c;
      res[i]->o=c+rcco*dco;
    } else {
      Vector nca=res[i]->ca-res[i]->n;
      Vector coca=res[i]->c-res[i]->ca;
      Vector coc=nca+coca;
      Vector ncoc=coc/coc.norm();
      Vector c=res[i]->c-ncoc*dco/2.0;
      Vector o=res[i]->c+ncoc*dco/2.0;
      res[i]->c=c;
      res[i]->o=o;
    }
   }

  // finish missing side chain atoms
  for (i=0; i<nres; i++) {
    if (!strcmp(res[i]->name(),"VAL")) {
      Val *sres=(Val *)res[i];
      sres->build(sres->cb,1.540,111.5,-120.0,sres->ca,sres->n,sres->c); 
      double dih=dihedral(sres->n,sres->ca,sres->cb,sres->cg1);
      sres->build(sres->cg1,1.521,110.5,dih-60.0,sres->cb,sres->ca,sres->n);
      sres->build(sres->cg2,1.521,110.5,dih+60.0,sres->cb,sres->ca,sres->n);
    }

    if (!strcmp(res[i]->name(),"LEU")) {
      Leu *sres=(Leu *)res[i];
      sres->build(sres->cb,1.530,110.5,-120.0,sres->ca,sres->n,sres->c); 
      Vector cbg=sres->cg-sres->cb;
      sres->cg=sres->cb+2.0*cbg;
      double dih=dihedral(sres->ca,sres->cb,sres->cg,sres->cd1);
      sres->build(sres->cd1,1.521,110.7,dih-60.0,sres->cg,sres->cb,sres->ca);
      sres->build(sres->cd2,1.521,110.7,dih+60.0,sres->cg,sres->cb,sres->ca);
    }
	
    if (!strcmp(res[i]->name(),"ILE")) {
      Ile *sres=(Ile *)res[i];
      sres->build(sres->cb,1.540,111.5,-120.0,sres->ca,sres->n,sres->c); 
      double dih=dihedral(sres->n,sres->ca,sres->cb,sres->cg1);
      sres->build(sres->cg2,1.530,110.4,dih-60.0,sres->cb,sres->ca,sres->n);
      sres->build(sres->cg1,1.521,110.5,dih+60.0,sres->cb,sres->ca,sres->n);
    }
	
    if (!strcmp(res[i]->name(),"PRO")) {
      Pro *sres=(Pro *)res[i];
      sres->build(sres->cb,1.530,104.0,-120.0,sres->ca,sres->n,sres->c); 
      Vector cgcd=(sres->cg*3.0-sres->cb)/2.0;
      if (i>0 && res[i-1]->pdbnum+1 == res[i]->pdbnum) {
        Vector cn=res[i-1]->c-sres->n;   
        Vector ncn=cn/cn.norm();
        Vector can=sres->ca-sres->n;
        Vector ncan=can/can.norm();
        
        double cosa=cos(125.0/180.0*3.141519265);
        double cosb=cos(110.0/180.0*3.141519265);
        double cosd=cos(125.0/180.0*3.141519265);
        
        double mu=cosa/(1-cosd)-cosb*cosd/(1-cosd);
        double lam=cosb-mu*cosd;
        Vector cdn=lam*ncn+mu*ncan;
        Vector ncdn=cdn/cdn.norm();
        sres->cd=sres->n+ncdn*1.46;
      } else {
        sres->build(sres->cd,1.46,110.0,-15.0,sres->n,sres->ca,sres->cb); 
      }
      sres->cg=sres->cd+(cgcd-sres->cd)*2.0;
    }
	
    if (!strcmp(res[i]->name(),"CYS")) {
      Cys *sres=(Cys *)res[i];
      sres->build(sres->cb,1.530,110.5,-120.0,sres->ca,sres->n,sres->c); 
      Vector cbsg=sres->sg-sres->cb;
      sres->sg=sres->cb+2.0*cbsg;
    }
	
    if (!strcmp(res[i]->name(),"MET")) {
      Met *sres=(Met *)res[i];
      sres->build(sres->cb,1.530,110.5,-120.0,sres->ca,sres->n,sres->c); 
      Vector cbg=sres->cg-sres->cb;
      sres->cg=sres->cb+2.0*cbg;
    }

    if (!strcmp(res[i]->name(),"PHE")) {
      Phe *sres=(Phe *)res[i];
      sres->build(sres->cb,1.530,110.5,-120.0,sres->ca,sres->n,sres->c); 
      Vector cbg=sres->cg-sres->cb;
      sres->cg=sres->cb+2.0*cbg;
      sres->build(sres->cd1,1.39,90.0,0.0,sres->ce1,sres->ce2,sres->cb);
      sres->build(sres->cd2,1.39,90.0,0.0,sres->ce2,sres->ce1,sres->cb);
      sres->build(sres->cz,1.39,120.0,0.0,sres->ce2,sres->cd2,sres->cb);
    }

    if (!strcmp(res[i]->name(),"TYR")) {
      Tyr *sres=(Tyr *)res[i];
      sres->build(sres->cb,1.530,110.5,-120.0,sres->ca,sres->n,sres->c); 
      Vector cbg=sres->cg-sres->cb;
      sres->cg=sres->cb+2.0*cbg;
      sres->build(sres->cd1,1.39,90.0,0.0,sres->ce1,sres->ce2,sres->cb);
      sres->build(sres->cd2,1.39,90.0,0.0,sres->ce2,sres->ce1,sres->cb);
      sres->build(sres->cz,1.39,120.0,0.0,sres->ce2,sres->cd2,sres->cb);
    }

    if (!strcmp(res[i]->name(),"TRP")) {
      Trp *sres=(Trp *)res[i];
      sres->build(sres->cb,1.530,110.5,-120.0,sres->ca,sres->n,sres->c); 
      Vector cbg=sres->cg-sres->cb;
      sres->cg=sres->cb+2.0*cbg;
      sres->build(sres->cd2,1.40,60.0,0.0,sres->ce3,sres->cz2,sres->ne1);
      sres->build(sres->ce2,1.345,26.0,0.0,sres->ne1,sres->cz2,sres->ce3);
      sres->build(sres->cd1,1.318,110.5,0.0,sres->ne1,sres->ce2,sres->cd2);
      sres->build(sres->cz3,1.40,119.5,0.0,sres->ce3,sres->cd2,sres->ce2);
      sres->build(sres->ch2,1.40,118.6,0.0,sres->cz2,sres->ce2,sres->cd2);
    }

    if (!strcmp(res[i]->name(),"SER")) {
      Ser *sres=(Ser *)res[i];
      sres->build(sres->cb,1.530,110.5,-120.0,sres->ca,sres->n,sres->c); 
      Vector cbog=sres->og-sres->cb;
      sres->og=sres->cb+2.0*cbog;
    }

    if (!strcmp(res[i]->name(),"THR")) {
      Thr *sres=(Thr *)res[i];
      sres->build(sres->cb,1.540,111.5,-120.0,sres->ca,sres->n,sres->c); 
      Vector cbog=sres->og1-sres->cb;
      sres->og1=sres->cb+2.0*cbog;
    }

    if (!strcmp(res[i]->name(),"ASN")) {
      Asn *sres=(Asn *)res[i];

      double dcac=1.516;
      double dcn=1.328;
      double dco=1.231;
      double dcac2=dcac*dcac;
      double dcn2=dcn*dcn;
      Vector rcaco=sres->od1-sres->cb;
      Vector rcan=sres->nd2-sres->cb;
      double dcaco=rcaco.norm();
      double dcan=rcan.norm();
      Vector ncaco=rcaco/dcaco;
      Vector ncan=rcan/dcan;
      double dcaco2=dcaco*dcaco;
      double dcan2=dcan*dcan;
      double a=(dcan2+dcac2-dcn2)/(2.0*dcan);
      double h2=dcac2-a*a;
      double h=sqrt(h2);
      double cosa=ncaco*ncan; 
      double sina=sqrt(1.0-cosa*cosa);
      double x=h/sina;   
      double ab=sqrt(x*x-h2);
      double aa=a-ab;
      Vector c=sres->cb+ncan*aa+ncaco*x;
      Vector rcco=sres->od1-c;
      rcco/=rcco.norm();
      sres->cg=c;
      sres->od1=c+rcco*dco;
    }

    if (!strcmp(res[i]->name(),"GLN")) {
      Gln *sres=(Gln *)res[i];
      sres->build(sres->cb,1.530,110.5,-120.0,sres->ca,sres->n,sres->c); 
      Vector cbg=sres->cg-sres->cb;
      sres->cg=sres->cb+2.0*cbg;

      double dcac=1.516;
      double dcn=1.328;
      double dco=1.231;
      double dcac2=dcac*dcac;
      double dcn2=dcn*dcn;
      Vector rcaco=sres->oe1-sres->cg;
      Vector rcan=sres->ne2-sres->cg;
      double dcaco=rcaco.norm();
      double dcan=rcan.norm();
      Vector ncaco=rcaco/dcaco;
      Vector ncan=rcan/dcan;
      double dcaco2=dcaco*dcaco;
      double dcan2=dcan*dcan;
      double a=(dcan2+dcac2-dcn2)/(2.0*dcan);
      double h2=dcac2-a*a;
      double h=sqrt(h2);
      double cosa=ncaco*ncan; 
      double sina=sqrt(1.0-cosa*cosa);
      double x=h/sina;   
      double ab=sqrt(x*x-h2);
      double aa=a-ab;
      Vector c=sres->cg+ncan*aa+ncaco*x;
      Vector rcco=sres->oe1-c;
      rcco/=rcco.norm();
      sres->cd=c;
      sres->oe1=c+rcco*dco;
    }

    if (!strcmp(res[i]->name(),"LYS")) {
      Lys *sres=(Lys *)res[i];
      sres->build(sres->cb,1.530,110.5,-120.0,sres->ca,sres->n,sres->c); 
      Vector cbg=sres->cg-sres->cb;
      sres->cg=sres->cb+2.0*cbg;
    }

    if (!strcmp(res[i]->name(),"ARG")) {
      Arg *sres=(Arg *)res[i];
      sres->build(sres->cb,1.530,110.5,-120.0,sres->ca,sres->n,sres->c); 
      Vector cbg=sres->cg-sres->cb;
      sres->cg=sres->cb+2.0*cbg;

      Vector necz=sres->nh1-sres->ne+sres->nh2-sres->ne;
      Vector nnecz=necz/necz.norm();
      sres->cz=sres->ne+nnecz*1.335;
    }

    if (!strcmp(res[i]->name(),"HIS") || !strcmp(res[i]->name(),"HSD") || 
	!strcmp(res[i]->name(),"HSE") || !strcmp(res[i]->name(),"HSP")) {
      His *sres=(His *)res[i];
      sres->build(sres->cb,1.530,110.5,-120.0,sres->ca,sres->n,sres->c); 
      Vector cbg=sres->cg-sres->cb;
      sres->cg=sres->cb+2.0*cbg;
      sres->build(sres->ce1,1.31,35.75,180.0,sres->ne2,sres->nd1,sres->cg);
      sres->build(sres->cd2,1.31,75.0,0.0,sres->ne2,sres->nd1,sres->cg);
    }

    if (!strcmp(res[i]->name(),"ASP")) {
      Asp *sres=(Asp *)res[i];
      Vector cbcg=sres->od1-sres->cb+sres->od2-sres->cb;
      Vector ncbcg=cbcg/cbcg.norm();
      sres->cg=sres->cb+ncbcg*1.52;
      Vector od1=(sres->od1-sres->cg)*2.0+sres->cg;
      Vector od2=(sres->od2-sres->cg)*2.0+sres->cg;
      sres->od1=od1;
      sres->od2=od2;
    }

    if (!strcmp(res[i]->name(),"GLU")) {
      Glu *sres=(Glu *)res[i];
      sres->build(sres->cb,1.530,110.5,-120.0,sres->ca,sres->n,sres->c); 
      Vector cbg=sres->cg-sres->cb;
      sres->cg=sres->cb+2.0*cbg;
      Vector cgcd=sres->oe1-sres->cg+sres->oe2-sres->cg;
      Vector ncgcd=cgcd/cgcd.norm();
      sres->cd=sres->cg+ncgcd*1.52;
      Vector oe1=(sres->oe1-sres->cd)*2.0+sres->cd;
      Vector oe2=(sres->oe2-sres->cd)*2.0+sres->cd;
      sres->oe1=oe1;
      sres->oe2=oe2;
    }
  }



  int atominx=1;
  for (int i=0; i<nres; i++) {
    res[i]->print(atominx);
  }

  delete mol;
  delete rstart;
  
  return 0;
}

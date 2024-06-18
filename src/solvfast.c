// fast version 2024, Michael Feig (feig@msu.edu)
// cubic only, assumes solute is centered

// *** C++ code ***
//
// requires:
//   vector.h
//   pdb.h
//   field.h
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "vector.h"
#include "pdb.h"
#include "field.h"

#define ANINT(n) (n>=0?floor(n+0.5):-floor(0.5-n)) 

void usage() {
 fprintf(stderr,"usage:   solvate [options] PDBfile\n");
 fprintf(stderr,"options: [-box PDBfile] [-boxwidth value]\n");
 fprintf(stderr,"         [-cutoff value]\n");
 fprintf(stderr,"         [-solvcut value]\n");
 fprintf(stderr,"         [-ioncut value]\n");
 fprintf(stderr,"         [-ions NAME:num[=NAME:num]]\n");
 fprintf(stderr,"         [-fixbox xmin xmax ymin ymax zmin zmax]\n");
 fprintf(stderr,"         [-periodic]\n");
 fprintf(stderr,"         [-tip3p]\n");
 fprintf(stderr,"         [-verbose]\n");
 exit(1);
}

int verbose=0;

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

int main(int argc, char **argv) {
  srand(time(0));

  char pdbfile[256];
  double cutoff=9.0;
  double boxwidth=18.662;
  char boxfile[256];

  double solvcut=2.10;

  double ioncut=4.0;

  int maxions=10;
  int nions=0;
  char **ionname=new char*[maxions];
  int *ioncnt=new int[maxions];

  char tip3pname[12];
  strcpy(tip3pname,"TIP3");

  int periodic=0;

  strcpy(boxfile,"water.pdb");
  if (argc<2) {
    usage();
  }

  double boxminx=0.0;
  double boxmaxx=0.0;
  double boxminy=0.0;
  double boxmaxy=0.0;
  double boxminz=0.0;
  double boxmaxz=0.0;
  int fixbox=0;
  int tip3p=0;

  int i;
  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i],"-box")) {
      strcpy(boxfile,argv[++i]);
    } else if (!strcmp(argv[i],"-boxwidth")) {
      boxwidth=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-cutoff")) {
      cutoff=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-solvcut")) {
      solvcut=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-ioncut")) {
      ioncut=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-fixbox")) {
      fixbox=1;
      boxminx=atof(argv[++i]);
      boxmaxx=atof(argv[++i]);
      boxminy=atof(argv[++i]);
      boxmaxy=atof(argv[++i]);
      boxminz=atof(argv[++i]);
      boxmaxz=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-tip3p")) {
      tip3p=1;
    } else if (!strcmp(argv[i],"-ions")) {
      Field farg(argv[++i],"=");
      for (int ifarg=0; ifarg<farg.number(); ifarg++) {
        Field fion(farg[ifarg],":");
        if (nions<maxions) {
          ionname[nions]=new char[12];
          strcpy(ionname[nions],fion[0]);
          ioncnt[nions]=atoi(fion[1]);
          nions++;
        }
      }
    } else if (!strcmp(argv[i],"-verbose")) {
      verbose=1;
    } else if (!strcmp(argv[i],"-periodic")) {
      periodic=1;
    } else if (!strcmp(argv[i],"-help")) {
      usage();
    } else {
      strcpy(pdbfile,argv[i]);
    }
  }
  
  double solvcutsq=solvcut*solvcut;
  double extsolvcut=solvcut+1.5;
  double extsolvcutsq=extsolvcut*extsolvcut;

  double watcutsq=1.7*1.7;

  double cutsq=cutoff*cutoff;
  
  PDBEntry *solute,*solvent;
  int maxsolute=10000000;
  int maxsolvent=10000000;
  int maxsolvbox=10000000;

  solute=new PDBEntry[maxsolute];
  solvent=new PDBEntry[maxsolvent];

  int *solvinx=new int[maxsolvbox]; 
  Vector *solvadd=new Vector[maxsolvbox];
  int *solvedge=new int[maxsolvbox];
  int *solvion=new int[maxsolvbox];
  int nedge=0;

  int nsolute,nsolvent,nsolvbox;
  int nsoluteres,nsolventres;

  nsolvbox=0;

  int *soluteStart=new int[5000000];
  int *solventStart=new int[5000000];

  readPDB(boxfile,solvent,nsolvent,solventStart,nsolventres);
  readPDB(pdbfile,solute,nsolute,soluteStart,nsoluteres);

  if (verbose) {
    fprintf(stderr,"read %d atoms, %d residues from %s\n",nsolvent,nsolventres,boxfile);
    fprintf(stderr,"read %d atoms, %d residues from %s\n",nsolute,nsoluteres,pdbfile);
  }
  
  Vector cofm;
  for (i=0; i<nsolute; i++) {
    cofm+=solute[i].coordinates();
  }
  cofm/=nsolute;
  if (verbose) {
    fprintf(stderr,"center of mass: %lf %lf %lf\n",cofm.x(),cofm.y(),cofm.z());
  }

  int atominx=1;
  int resinx=solute[0].residueNumber()-1;

  int ix,iy,iz,in,ip,ipa,iwa;

  Vector min(99999.0,99999.0,99999.0), max(-99999.0,-99999.0,-99999.0);
  
  for (i=0; i<nsoluteres; i++) {
    for (ipa=soluteStart[i]; ipa<soluteStart[i+1]; ipa++) {
      PDBEntry ts(atominx++,solute[ipa].atomName(),solute[ipa].residueName(),
		  solute[ipa].residueNumber(),solute[ipa].coordinates(),
                  solute[ipa].segmentID());
      resinx=solute[ipa].residueNumber();
      ts.chain()=solute[ipa].chain();
      ts.write(stdout);
      if (solute[ipa].coordinates().x()<min.x()) min.x()=solute[ipa].coordinates().x();
      if (solute[ipa].coordinates().x()>max.x()) max.x()=solute[ipa].coordinates().x();
      if (solute[ipa].coordinates().y()<min.y()) min.y()=solute[ipa].coordinates().y();
      if (solute[ipa].coordinates().y()>max.y()) max.y()=solute[ipa].coordinates().y();
      if (solute[ipa].coordinates().z()<min.z()) min.z()=solute[ipa].coordinates().z();
      if (solute[ipa].coordinates().z()>max.z()) max.z()=solute[ipa].coordinates().z();
    }
  }
  printf("TER\n");

  if (fixbox) {
    min=Vector(boxminx,boxminy,boxminz);
    max=Vector(boxmaxx,boxmaxy,boxmaxz);
  } else {
    if (cofm.x()-min.x()>max.x()-cofm.x()) 
      max.x()=cofm.x()+cofm.x()-min.x();
    else 
      min.x()=cofm.x()-(max.x()-cofm.x());

    if (cofm.y()-min.y()>max.y()-cofm.y()) 
      max.y()=cofm.y()+cofm.y()-min.y();
    else 
      min.y()=cofm.y()-(max.y()-cofm.y());

    if (cofm.z()-min.z()>max.z()-cofm.z()) 
      max.z()=cofm.z()+cofm.z()-min.z();
    else 
      min.z()=cofm.z()-(max.z()-cofm.z());
    
    Vector cVec;
    if (max.x()>=max.y() && max.x()>=max.z()) 
      cVec=Vector(max.x(),max.x(),max.x());
    else if (max.y()>=max.x() && max.y()>=max.z()) 
      cVec=Vector(max.y(),max.y(),max.y());
    else 
      cVec=Vector(max.z(),max.z(),max.z());
    max=cVec;

    if (min.x()<=min.y() && min.x()<=min.z()) 
      cVec=Vector(min.x(),min.x(),min.x());
    else if (min.y()<=min.x() && min.y()<=min.z()) 
      cVec=Vector(min.y(),min.y(),min.y());
    else 
      cVec=Vector(min.z(),min.z(),min.z());
    min=cVec;

    min-=cutoff;
    max+=cutoff;
  }

  Vector dim=max-min;

  if (verbose) {
    fprintf(stderr,"min: %lf %lf %lf\n",min.x(),min.y(),min.z());
    fprintf(stderr,"max: %lf %lf %lf\n",max.x(),max.y(),max.z());
    fprintf(stderr,"dim: %lf %lf %lf\n",dim.x(),dim.y(),dim.z());
  }

  double dxgrid=solvcut;
  int ngridx=int(dim.x()/dxgrid)+1;
  int ngridy=int(dim.y()/dxgrid)+1;
  int ngridz=int(dim.z()/dxgrid)+1;

  fprintf(stderr,"lookup grid: %d x %d x %d\n",ngridx,ngridy,ngridz);

  int nmax=100;
  int *lookup=new int[ngridx*ngridy*ngridz*nmax];
  int *nlookup=new int[ngridx*ngridy*ngridz];
  for (int i=0; i<ngridx*ngridy*ngridz; i++) {
    nlookup[i]=0;
  }
  for (ip=0; ip<nsolute; ip++) {
    Vector xyz=solute[ip].coordinates();
    if (periodic) {
      xyz-=Vector(ANINT(xyz.x()/dim.x())*dim.x(),ANINT(xyz.y()/dim.y())*dim.y(),ANINT(xyz.z()/dim.z())*dim.z());
    }
//    Vector d=solute[ip].coordinates()-min;
    Vector d=xyz-min;
    int igx=int(d.x()/dxgrid);  
    int igy=int(d.y()/dxgrid);  
    int igz=int(d.z()/dxgrid);  
    for (ix=igx-1; ix<=igx+1; ix++) {
      int iix=ix;
      if (periodic && iix<0) iix+=ngridx;
      if (periodic && iix>=ngridx) iix-=ngridx;
      if (iix>=0 && iix<ngridx) {
        for (iy=igy-1; iy<=igy+1; iy++) {
          int iiy=iy;
          if (periodic && iiy<0) iiy+=ngridy;
          if (periodic && iiy>=ngridy) iiy-=ngridy;
          if (iiy>=0 && iiy<ngridy) {
            for (iz=igz-1; iz<=igz+1; iz++) {
              int iiz=iz;
              if (periodic && iiz<0) iiz+=ngridz;
              if (periodic && iiz>=ngridz) iiz-=ngridz;
              if (iiz>=0 && iiz<ngridz) {
                int inx=iiz*ngridx*ngridy+iiy*ngridx+iix;
                if (nlookup[inx]<nmax) {
                 lookup[inx+nlookup[inx]*ngridx*ngridy*ngridz]=ip;
                 nlookup[inx]++;
                } else {
                 fprintf(stderr,"reached limit %d %d %d\n",iix,iiy,iiz);
                 exit(1);
                }
              }
            }
          }
        }
      }
    }
  } 

  int xmult=int(dim.x()/boxwidth)*2;
  int ymult=int(dim.y()/boxwidth)*2;
  int zmult=int(dim.z()/boxwidth)*2;
 
  fprintf(stderr,"box size: %lf x %lf x %lf\n",dim.x(),dim.y(),dim.z());

  for (ix=0; ix<=xmult; ix++) {
    for (iy=0; iy<=ymult; iy++) {
      for (iz=0; iz<=zmult; iz++) {
	Vector a=Vector((double)ix,(double)iy,(double)iz);
	a*=boxwidth;
        a+=min; 
	for (in=0; in<nsolventres; in++) {
	  int good=0;
	  Vector c=solvent[solventStart[in]].coordinates()+a;
	  if (c.x()>min.x() && c.x()<max.x() &&
	      c.y()>min.y() && c.y()<max.y() &&
	      c.z()>min.z() && c.z()<max.z()) {
	    good=1;

            int edge=0;
            if (c.x()<min.x()+2.0 || c.x()>max.x()-2.0 || 
                c.y()<min.y()+2.0 || c.y()>max.y()-2.0 ||
                c.z()<min.z()+2.0 || c.z()>max.z()-2.0) {
              edge=1;
            } 

 	    for (iwa=solventStart[in]; iwa<solventStart[in+1] && good; iwa++) {
              Vector cc=solvent[iwa].coordinates()+a;
              Vector ta=cc-min;
              int igx=int(ta.x()/dxgrid);  
              int igy=int(ta.y()/dxgrid);  
              int igz=int(ta.z()/dxgrid);  
              int inx=igz*ngridx*ngridy+igy*ngridx+igx;

              if (inx>=0 && inx<ngridx*ngridy*ngridz && nlookup[inx]>0) {
  	        for (ip=0; ip<nlookup[inx] && good; ip++) {
                  Vector svec=solute[lookup[inx+ip*ngridx*ngridy*ngridz]].coordinates();
	          Vector d=cc-svec;
                  if (periodic) {
                     d-=Vector(ANINT(d.x()/dim.x())*dim.x(),ANINT(d.y()/dim.y())*dim.y(),ANINT(d.z()/dim.z())*dim.z());
                  }
	  	  double dval=d*d;
		  if (dval<solvcutsq) {
		    good=0;
                  }
                }
              }
            }

	    if (good && edge) {
              for (int iedge=0; iedge<nedge && good; iedge++) {
                int inxedge=solvinx[solvedge[iedge]];
                Vector sadd=solvadd[solvedge[iedge]];
                for (int iwah=solventStart[inxedge]; iwah<solventStart[inxedge+1] && good; iwah++) {
                  for (int iwan=solventStart[in]; iwan<solventStart[in+1] && good; iwan++) {
                    Vector dww=(solvent[iwan].coordinates()+a)-(solvent[iwah].coordinates()+sadd);
                    dww.pbc(dim);
                    double dwwval=dww*dww;
                    if (dwwval<watcutsq) good=0;
                  }
                }
              }
            }
              
            if (good) { 
              solvinx[nsolvbox]=in;
              solvion[nsolvbox]=0;
              solvadd[nsolvbox]=a;
              if (edge) {
                if (nedge>maxsolvbox) { 
                   fprintf(stderr,"exceeded maximum number of edge solvents (%d)\n",maxsolvbox);
                   exit(1);
                }
                solvedge[nedge++]=nsolvbox; 
              }
              nsolvbox++;
              if (verbose && nsolvbox%10000==0) { 
                fprintf(stderr,"added %d waters\n",nsolvbox);
              }
	    }
	  }
	}
      }
    }
  }

  int nseg=0;
  char segname[12];

  if (nions>0) {
    if (verbose) {
      fprintf(stderr,"adding ions\n");
    }
    for (int ii=0; ii<nions; ii++) {
      if (verbose) {
        fprintf(stderr,"%d ions type >%s<\n",ioncnt[ii],ionname[ii]);
      }
      resinx=1;
      for (int inum=0; inum<ioncnt[ii]; inum++) {
        if (resinx==1) sprintf(segname,"I%03d",nseg);
        int irnd;
        int inx;
        do {
          irnd=rand()%nsolvbox;
          Vector ta=solvent[solventStart[solvinx[irnd]]].coordinates()+solvadd[irnd]-min;
          int igx=int(ta.x()/dxgrid);  
          int igy=int(ta.y()/dxgrid);  
          int igz=int(ta.z()/dxgrid);  
          inx=igz*ngridx*ngridy+igy*ngridx+igx;
        } while (solvion[irnd] || nlookup[inx]>0);
        solvion[irnd]=1; 
        PDBEntry t(atominx,ionname[ii],ionname[ii],
                   resinx,solvent[solventStart[solvinx[irnd]]].coordinates()+solvadd[irnd],
                   segname);
        t.write(stdout);
        if (atominx<9999999) atominx++;
        resinx++;
        if (resinx>=10000) {
          printf("TER\n");
          resinx=1;
          nseg++;
        }
      }
      printf("TER\n");
      nseg++; 
    }
  }

  if (verbose) {
    fprintf(stderr,"writing out waters\n");
  }
  resinx=1;
  nseg=0;
  for (int isb=0; isb<nsolvbox; isb++) {
    if (!solvion[isb]) {
      if (resinx==1) sprintf(segname,"W%03d",nseg);

      int in=solvinx[isb];
      Vector a=solvadd[isb];
      for (iwa=solventStart[in]; iwa<solventStart[in+1]; iwa++) {
         PDBEntry t(atominx,solvent[iwa].atomName(),(tip3p?tip3pname:solvent[iwa].residueName()),
                    resinx,solvent[iwa].coordinates()+a,segname);
         t.write(stdout);
         if (atominx<9999999) atominx++;
      }
      resinx++;
      if (resinx>=10000) {
        printf("TER\n");
        resinx=1;
        nseg++;
      }
    }
  } 

  if (resinx>1) printf("TER\n");
  printf("END\n");

  delete solvinx;
  delete solvadd;
  delete solute;
  delete solvent;
  delete soluteStart;
  delete solventStart;
}

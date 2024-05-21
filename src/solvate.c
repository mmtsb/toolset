// http://mmtsb.scripps.edu/doc/solvate.html
// 2000, Michael Feig (meikel@scripps.edu)
// Brooks group, The Scripps Research Institute 
// rewritten 2023, Michael Feig (feig@msu.edu)

// *** C++ code ***
//
// requires:
//   vector.h
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "vector.h"
#include "pdb.h"

void usage() {
 fprintf(stderr,"usage:   solvate [options] PDBfile\n");
 fprintf(stderr,"options: [-box PDBfile] [-boxwidth value]\n");
 fprintf(stderr,"         [-cutoff value]\n");
 fprintf(stderr,"         [-res min max] [-cubic] [-octahedron]\n");
 fprintf(stderr,"         [-solvcut value]\n");
 fprintf(stderr,"         [-[no]center]\n");
 fprintf(stderr,"         [-fixbox xmin xmax ymin ymax zmin zmax]\n");
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
  char pdbfile[256];
  double cutoff=9.0;
  double boxwidth=18.662;
  char boxfile[256];
  int resmin=-1;
  int resmax=-1;

  int cubicmode=0;
  int octamode=0;
  int center=1;

  double solvcut=2.10;

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

  int i;
  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i],"-cubic")) {
      cubicmode=1;
    } else if (!strcmp(argv[i],"-octahedron")) {
      octamode=1;
      cubicmode=1;
    } else if (!strcmp(argv[i],"-box")) {
      strcpy(boxfile,argv[++i]);
    } else if (!strcmp(argv[i],"-res")) {
      resmin=atoi(argv[++i])-1;
      resmax=atoi(argv[++i])-1;
    } else if (!strcmp(argv[i],"-boxwidth")) {
      boxwidth=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-cutoff")) {
      cutoff=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-solvcut")) {
      solvcut=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-center")) {
      center=1;
    } else if (!strcmp(argv[i],"-nocenter")) {
      center=0;
    } else if (!strcmp(argv[i],"-fixbox")) {
      fixbox=1;
      boxminx=atof(argv[++i]);
      boxmaxx=atof(argv[++i]);
      boxminy=atof(argv[++i]);
      boxmaxy=atof(argv[++i]);
      boxminz=atof(argv[++i]);
      boxmaxz=atof(argv[++i]);
      center=0;
      octamode=0;
    } else if (!strcmp(argv[i],"-verbose")) {
      verbose=1;
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
      if (center) {
        solute[ipa].coordinates()-=cofm;
      }
      PDBEntry ts(atominx++,solute[ipa].atomName(),solute[ipa].residueName(),
		  solute[ipa].residueNumber(),solute[ipa].coordinates(),
                  solute[ipa].segmentID());
      resinx=solute[ipa].residueNumber();
      ts.chain()=solute[ipa].chain();
      ts.write(stdout);

      if (solute[ipa].coordinates().x()<min.x()) {
        min.x()=solute[ipa].coordinates().x();
      }
      if (solute[ipa].coordinates().x()>max.x()) {
        max.x()=solute[ipa].coordinates().x();
      }
      if (solute[ipa].coordinates().y()<min.y()) {
	min.y()=solute[ipa].coordinates().y();
      }
      if (solute[ipa].coordinates().y()>max.y()) {
	max.y()=solute[ipa].coordinates().y();
      }
      if (solute[ipa].coordinates().z()<min.z()) {
	min.z()=solute[ipa].coordinates().z();
      }
      if (solute[ipa].coordinates().z()>max.z()) {
	max.z()=solute[ipa].coordinates().z();
      }
    }
  }

  if (fixbox) {
    min=Vector(boxminx,boxminy,boxminz);
    max=Vector(boxmaxx,boxmaxy,boxmaxz);
  } else {
    if (cofm.x()-min.x()>max.x()-cofm.x()) {
      max.x()=cofm.x()+cofm.x()-min.x();
    } else {
      min.x()=cofm.x()-(max.x()-cofm.x());
    }

    if (cofm.y()-min.y()>max.y()-cofm.y()) {
      max.y()=cofm.y()+cofm.y()-min.y();
    } else {
      min.y()=cofm.y()-(max.y()-cofm.y());
    }

    if (cofm.z()-min.z()>max.z()-cofm.z()) {
      max.z()=cofm.z()+cofm.z()-min.z();
    } else {
      min.z()=cofm.z()-(max.z()-cofm.z());
    }
    min-=cutoff;
    max+=cutoff;

    if (cubicmode) {
      Vector cVec;
      if (max.x()>=max.y() && max.x()>=max.z()) {
        cVec=Vector(max.x(),max.x(),max.x());
      } else if (max.y()>=max.x() && max.y()>=max.z()) {
        cVec=Vector(max.y(),max.y(),max.y());
      } else {
        cVec=Vector(max.z(),max.z(),max.z());
      }
      max=cVec;
      if (min.x()<=min.y() && min.x()<=min.z()) {
        cVec=Vector(min.x(),min.x(),min.x());
      } else if (min.y()<=min.x() && min.y()<=min.z()) {
        cVec=Vector(min.y(),min.y(),min.y());
      } else {
        cVec=Vector(min.z(),min.z(),min.z());
      }
      min=cVec;
    }
  }

  Vector dim=max-min;

  if (verbose) {
    fprintf(stderr,"min: %lf %lf %lf\n",min.x(),min.y(),min.z());
    fprintf(stderr,"max: %lf %lf %lf\n",max.x(),max.y(),max.z());
    fprintf(stderr,"dim: %lf %lf %lf\n",dim.x(),dim.y(),dim.z());
  }

  int xmult=int(dim.x()/boxwidth)*2;
  int ymult=int(dim.y()/boxwidth)*2;
  int zmult=int(dim.z()/boxwidth)*2;
 
  if (octamode) {
    fprintf(stderr,"box size: %lf x %lf x %lf\n",
            dim.x()*sqrt(3.0/4.0),dim.y()*sqrt(3.0/4.0),dim.z()*sqrt(3.0/4.0));
    fprintf(stderr,"surrounding cube: %lf x %lf x %lf\n",dim.x(),dim.y(),dim.z());
  } else {
    fprintf(stderr,"box size: %lf x %lf x %lf\n",dim.x(),dim.y(),dim.z());
  } 

  for (ix=0; ix<=xmult; ix++) {
    for (iy=0; iy<=ymult; iy++) {
      for (iz=0; iz<=zmult; iz++) {
	Vector a=Vector((double)ix,(double)iy,(double)iz);
	a*=boxwidth;
        a+=min; 
	for (in=0; in<nsolventres; in++) {
	  int good=0;
	  int goodcut=0;
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
            //fprintf(stderr,"have water at %lf %lf %lf\n",c.x(),c.y(),c.z());

	    if (resmin>=0 && resmax>=resmin) {
	      for (ip=0; ip<nsoluteres && good; ip++) {
		for (ipa=soluteStart[ip]; ipa<soluteStart[ip+1] && good; ipa++) {
		  Vector d=c-solute[ipa].coordinates();
                  d.pbc(dim);
		  double dval=d*d;
		  if (dval<solvcutsq) {
		    good=0;
		  } else if (dval<extsolvcutsq) {
		    for (iwa=solventStart[in]+1; iwa<solventStart[in+1] && good; iwa++) {
		      Vector dw=solvent[iwa].coordinates()+a-solute[ipa].coordinates();
                      dw.pbc(dim);
		      double dwval=dw*dw;
		      
		      if (dwval<solvcutsq) {
			good=0;
		      } else if (ip>=resmin && ip<=resmax && dwval<cutsq) {
			goodcut++;
		      }
		    }
		  } else if (ip>=resmin && ip<=resmax && dval<cutsq) {
		    goodcut++;
		  }
		}
	      }
	    } else {
	      goodcut=1;
	      for (ip=0; ip<nsoluteres && good; ip++) {
		for (ipa=soluteStart[ip]; ipa<soluteStart[ip+1] && good; ipa++) {
		  Vector d=c-solute[ipa].coordinates();
                  d.pbc(dim);
		  double dval=d*d;
		  if (dval<solvcutsq) {
		    good=0;
		  } else if (dval<extsolvcutsq) {
		    for (iwa=solventStart[in]+1; iwa<solventStart[in+1] && good; iwa++) {
		      Vector dw=solvent[iwa].coordinates()+a-solute[ipa].coordinates();
                      dw.pbc(dim);
		      double dwval=dw*dw;
		      if (dwval<solvcutsq) {
			good=0;
		      }
		    }
		  }
		}
	      }
	    }

	    iwa=solventStart[in];
            double tx=solvent[iwa].coordinates().x()+a.x();
	    double ty=solvent[iwa].coordinates().y()+a.y();
            double tz=solvent[iwa].coordinates().z()+a.z();
	    double rcf=(fabs(tx)+fabs(ty)+fabs(tz));
	    if (good && goodcut>0 && (!octamode || rcf<1.5*max.x())) {
              if (edge) {
                for (int ihave=0; ihave<nsolvbox && good; ihave++) {
                  if (solvedge[ihave]) {
                    for (int iwah=solventStart[solvinx[ihave]]; iwah<solventStart[solvinx[ihave]+1] && good; iwah++) {
                      for (int iwan=solventStart[in]; iwan<solventStart[in+1] && good; iwan++) {
                        Vector dww=(solvent[iwan].coordinates()+a)-(solvent[iwah].coordinates()+solvadd[ihave]);
                        dww.pbc(dim);
                        double dwwval=dww*dww;
                        if (dwwval<watcutsq) {
                          good=0;
                        }
                      }
                    }
                  }
                }
              }

              if (good) { 
                solvinx[nsolvbox]=in;
                solvedge[nsolvbox]=edge;
                solvadd[nsolvbox++]=a;
                if (verbose && nsolvbox%1000==0) { 
                  fprintf(stderr,"added %d waters\n",nsolvbox);
                }
              }
	    }
	  }
	}
      }
    }
  }

  for (int isb=0; isb<nsolvbox; isb++) {
      resinx++;
      int in=solvinx[isb];
      Vector a=solvadd[isb];
      for (iwa=solventStart[in]; iwa<solventStart[in+1]; iwa++) {
         PDBEntry t(atominx++,solvent[iwa].atomName(),solvent[iwa].residueName(),
                    resinx,solvent[iwa].coordinates()+a);
         t.write(stdout);
      }
  } 

  printf("END\n");

  delete solvinx;
  delete solvadd;
  delete solute;
  delete solvent;
  delete soluteStart;
  delete solventStart;
}

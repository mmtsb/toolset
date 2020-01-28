// http://mmtsb.scripps.edu/doc/kclust.html
// 2000, Michael Feig, Brooks group, TSRI
//

// *** C++ code ***
//
// requires:
//   vector.h pdb.h pdb.c field.h field.c
//   charmmlsq.f
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "field.h"
#include "vector.h"
#include "pdb.h"

void usage() {
  fprintf(stderr,"usage:   kclust [options] [inputFile]\n");
  fprintf(stderr,"options: [-pdb | -sicho]\n");
  fprintf(stderr,"         [-mode rmsd|phipsi|phi|psi|mix mixfactor]\n");
  fprintf(stderr,"         [-radius value]\n");
  fprintf(stderr,"         [-maxerr value] [-iterate]\n");
  fprintf(stderr,"         [-centroid] [-cdist]\n");
  fprintf(stderr,"         [-ca | -cb | -cab | -heavy | -all]\n");
  fprintf(stderr,"         [-l min:max[=min:max=...]]\n");
  fprintf(stderr,"         [-lsqfit] [-fitxl] [-fit min:max[=min:max=...]]\n");
  exit(1);
}

enum InputModeEnum { PDBMODE, SICHOMODE };
enum ClusterModeEnum { RMSDMODE, PHIPSIMODE, PHIMODE, PSIMODE, MIXMODE };

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

const int MAXATOMS=40000;
const int MAXRES=4000;

const int MAX_S=50000;

int verbose=0;

int dofit[MAXRES];
int dormsd[MAXRES];

PDBEntry pdbRef[MAXATOMS];

int lsqfitFlag=0;

double mixfactor=0.5;

ClusterModeEnum clusterMode=RMSDMODE;

// *** Structure ***********************************

class Structure {
 private:
  Vector *crd;
  double *phi;
  double *psi;

  int *alt;

  int nres;
  int natoms;
  int inx;

  int *resLookup;

  int findAlternate(PDBEntry *pdbs, int i, int natoms, char *aname);

 public:
  Structure(int i=0);
  Structure(const Structure& conf);
  ~Structure();
  
  Structure& operator=(const Structure& conf);

  Vector& coordinates(int inx) { return crd[inx]; }
  double& phiDihedral(int inx) { return phi[inx]; }
  double& psiDihedral(int inx) { return psi[inx]; }

  int residues() const { return nres; }
  int atoms() const { return natoms; }
  int index() const { return inx; }
  
  void read(char *fname, InputModeEnum=PDBMODE, 
	    SelEnum sel=ALL, PDBEntry *pdb=0);

  void calcDihedrals();

  double rmsdcoor(Structure& s, int lsqFlag);
  double rmsdphi(Structure& s);
  double rmsdpsi(Structure& s);

  double rmsd(Structure& s);

  void lsqFit(Structure& s);

  void add(Structure &s, int n, int lsqfit);

  double dihedral(int i, int j, int k, int l);
};

Structure::Structure(int i) : 
  natoms(0), nres(0), inx(i), crd(0), phi(0), psi(0), resLookup(0), alt(0) {}

Structure::~Structure() {
  if (crd!=0) 
    delete[] crd;

  if (alt!=0) 
    delete[] alt;

  if (resLookup!=0)
    delete resLookup;

  if (phi!=0)
    delete phi;

  if (psi!=0)
    delete psi;
}

Structure::Structure(const Structure& conf) :
  natoms(conf.natoms), nres(conf.nres), inx(conf.inx) {

  int i;
  if (conf.crd!=0) {
    crd=new Vector[natoms];
    for (i=0; i<natoms; i++)
      crd[i]=conf.crd[i];
  } else 
    crd=0;

  if (conf.alt!=0) {
    alt=new int[natoms];
    for (i=0; i<natoms; i++) 
      alt[i]=conf.alt[i];
  } else 
    alt=0;
  
  if (conf.phi!=0) {
    phi=new double[nres];
    for (i=0; i<nres; i++)
      phi[i]=conf.phi[i];
  } else
    phi=0;

  if (conf.psi!=0) {
    psi=new double[nres];
    for (i=0; i<nres; i++)
      psi[i]=conf.psi[i];
  } else
    psi=0;

  if (conf.resLookup!=0) {
    resLookup=new int[nres];
    for (i=0; i<nres; i++)
      resLookup[i]=conf.resLookup[i];
  } else 
    resLookup=0;
}

Structure& Structure::operator=(const Structure& conf) {
  natoms=conf.natoms;
  nres=conf.nres;
  inx=conf.inx;

  if (crd!=0)
    delete[] crd;

  if (alt!=0)
    delete[] alt;

  int i;

  crd=new Vector[natoms];
  for (i=0; i<natoms; i++)
    crd[i]=conf.crd[i];

  alt=new int[natoms];
  for (i=0; i<natoms; i++)
    alt[i]=conf.alt[i];

  if (conf.phi!=0) {
    if (phi!=0)
      delete[] phi;
    phi=new double[nres];
    for (i=0; i<nres; i++)
      phi[i]=conf.phi[i];
  }

  if (conf.psi!=0) {
    if (psi!=0) 
      delete[] psi;
    psi=new double[nres];
    for (i=0; i<nres; i++)
      psi[i]=conf.psi[i];
  }

  if (conf.resLookup!=0) {
    if (resLookup!=0)
      delete[] resLookup;
    resLookup=new int[nres];
    for (i=0; i<nres; i++)
      resLookup[i]=conf.resLookup[i];
  } 

  return (*this);
}

void Structure::add(Structure &s, int n, int lsqfit) {
  int i;

  if (lsqfit) s.lsqFit(*this);

  for (i=0; i<natoms; i++)
    crd[i]=((double)n*crd[i]+s.crd[i])/(n+1);
  
  if (phi!=0)
    for (i=0; i<nres; i++)
      phi[i]=((double)n*phi[i]+s.phi[i])/(double)(n+1);

  if (psi!=0) 
    for (i=0; i<nres; i++)
      psi[i]=((double)n*psi[i]+s.psi[i])/(double)(n+1);
}

int Structure::findAlternate(PDBEntry *pdbs, int i, int natoms, char *aname) {
  int resnum=pdbs[i].residueNumber();

  for (int j=i-1; j>=0 && pdbs[j].residueNumber()==resnum; j--) {
    if (!strcmp(pdbs[j].atomName(),aname)) {
      return j;
    }
  }

  for (int j=i+1; j<natoms && pdbs[j].residueNumber()==resnum; j++) {
    if (!strcmp(pdbs[j].atomName(),aname)) {
      return j;
    }
  }

  return -1;
}


void Structure::read(char *fname, InputModeEnum mode, 
		     SelEnum sel, PDBEntry *pdbs) {
  static int lastN=0;

  Vector *v=new Vector[MAXATOMS];
  
  FILE *fptr=0;
  int closepipe=0;

  if (!access(fname,R_OK)) {
    fptr=fopen(fname,"r");
  } else {
    char zipfname[1024];
    sprintf(zipfname,"%s.gz",fname);
    if (!access(zipfname,R_OK)) {
      char cmd[1024];
      sprintf(cmd,"gunzip -c %s",zipfname);
      fptr=popen(cmd,"r");
      closepipe=1;
    }
  }  

  if (fptr==0) {
    fprintf(stderr,"cannot open file %s\n",fname);
    exit(1);
  }

  PDBEntry pdb;

  int num=-1;
  int n=0;
  while (!feof(fptr)) {
    if (mode==PDBMODE) {
      if (pdb.read(fptr,sel)>0) {
	if (pdbs[n].atomIndex()<=0)
	  pdbs[n]=pdb;
	v[n++]=pdb.coordinates();
      }
    } else if (mode==SICHOMODE) {
      if (n==0) {
	fscanf(fptr,"%d",&num);
	n++;
      } else {
	int ix,iy,iz;
	if (fscanf(fptr,"%5d%5d%5d",&ix,&iy,&iz)>0) {
          if (n>1 && n<num) 
	    v[n-2]=Vector((double)ix,(double)iy,(double)iz);
	  n++;
	} 
      }
    }
  }

  if (closepipe)
    pclose(fptr);
  else
    fclose(fptr);

  if (mode==SICHOMODE) 
    n-=3;

  if (lastN==0) 
    lastN=n;
  else if (lastN!=n) {
    fprintf(stderr,
	    "number of atoms do not match in file %s (%d, expected: %d)\n",
	    fname,n,lastN);
    exit(1);
  }

  if (natoms==0) {
    natoms=n;
    crd=new Vector[natoms];
    phi=new double[natoms];
    psi=new double[natoms];
    alt=new int[natoms];
  } else if (natoms!=n) {
    fprintf(stderr,"number of atoms do not match (file: %d, init: %d)\n",n,natoms);
    exit(1);
  }
  
  for (int i=0; i<natoms; i++) {
    crd[i]=v[i];

    if (!strcmp(pdbs[i].atomName(),"C7") && *pdbs[i].residueName()=='W') {
      alt[i]=findAlternate(pdbs,i,natoms,"C12");
    } else if (!strcmp(pdbs[i].atomName(),"C12") && *pdbs[i].residueName()=='W') {
      alt[i]=findAlternate(pdbs,i,natoms,"C7");
    } else if (!strcmp(pdbs[i].atomName(),"C10") && *pdbs[i].residueName()=='W') {
      alt[i]=findAlternate(pdbs,i,natoms,"C8");
    } else if (!strcmp(pdbs[i].atomName(),"C8") && *pdbs[i].residueName()=='W') {
      alt[i]=findAlternate(pdbs,i,natoms,"C10");
    } else {
      alt[i]=-1;
    }
  }

  delete[] v;
}

void Structure::lsqFit(Structure& s) {
  Vector tc1(0,0,0);
  Vector tc2(0,0,0);
    
  int i;
  int nlsqfit=0;
  for (i=0; i<natoms; i++) {
    if (dofit[pdbRef[i].residueNumber()] &&
	pdbRef[i].selected()) {
      tc1+=crd[i];
      tc2+=s.crd[i];
      nlsqfit++;
    }
  }
  tc1/=(double)nlsqfit;
  tc2/=(double)nlsqfit;

  double r[9];
  double u[9];
  for (i=0; i<9; i++) 
    r[i]=0.0;
  
  for (i=0; i<natoms; i++) {
    crd[i]-=tc1;
    if (dofit[pdbRef[i].residueNumber()] &&
	pdbRef[i].selected()) {
      Vector ref=s.crd[i]-tc2;

      r[0]+=crd[i].x()*ref.x();
      r[1]+=crd[i].x()*ref.y();
      r[2]+=crd[i].x()*ref.z();
      
      r[3]+=crd[i].y()*ref.x();
      r[4]+=crd[i].y()*ref.y();
      r[5]+=crd[i].y()*ref.z();
      
      r[6]+=crd[i].z()*ref.x();
      r[7]+=crd[i].z()*ref.y();
      r[8]+=crd[i].z()*ref.z();
    }
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
  
  for (i=0; i<natoms; i++) {
    Vector nv(u[0]*crd[i].x()+u[3]*crd[i].y()+u[6]*crd[i].z(),
	      u[1]*crd[i].x()+u[4]*crd[i].y()+u[7]*crd[i].z(),
	      u[2]*crd[i].x()+u[5]*crd[i].y()+u[8]*crd[i].z());
    crd[i]=nv+tc2;
  }
}

double Structure::rmsdphi(Structure& s) {
  register int i;
  register double d=0.0;
  register int nd=0;

  //  fprintf(stderr,"nres: %d\n",nres);
  for (i=0; i<nres; i++) {
    if (dormsd[resLookup[i]]) {
      double dv=phi[i]-s.phi[i];
      d+=dv*dv;
      nd++;
    }
  }
  //  fprintf(stderr,"phi RMSD: nd: %d d: %lf, %lf\n",nd,d,sqrt(d/(double)nd));
  return sqrt(d/(double)nd);
}

double Structure::rmsdpsi(Structure& s) {
  register int i;
  register double d=0.0;
  register int nd=0;

  for (i=0; i<nres; i++) {
    if (dormsd[resLookup[i]]) {
      double dv=psi[i]-s.psi[i];
      d+=dv*dv;
      nd++;
    }
  }
  return sqrt(d/(double)nd);
}

double Structure::rmsdcoor(Structure& s, int lsqFlag) {
  register int i;
  register double d=0.0;
  register int nd=0;

  if (lsqFlag) {
    Structure ns(s);
    ns.lsqFit(*this);
    for (i=0; i<natoms; i++) {
      //      fprintf(stderr," atom %d, res: %d dormsd: %d, selected: %d\n",
      //	      i,pdbRef[i].residueNumber(),dormsd[pdbRef[i].residueNumber()],
      //	      pdbRef[i].selected());
      if (dormsd[pdbRef[i].residueNumber()] &&
	  pdbRef[i].selected()) {
	Vector dv;
	dv=crd[i]-ns.crd[i];

	if (alt[i]>=0) {
	  double d1=dv*dv;
	  dv=crd[i]-ns.crd[alt[i]];
	  double d2=dv*dv;

	  if (d1<d2) {
	    d+=d1;
	  } else {
	    d+=d2;
	  }
	} else {
	  d+=dv*dv;
	}

	nd++;
	//	fprintf(stderr," atom %d, res: %d, this: %lf %lf %lf, cmp: %lf %lf %lf\n",
	//		i,pdbRef[i].residueNumber(),crd[i].x(),crd[i].y(),crd[i].z(),
	//		ns.crd[i].x(),ns.crd[i].y(),ns.crd[i].z());
      }
    }
  } else {
    for (i=0; i<natoms; i++) {
      if (dormsd[pdbRef[i].residueNumber()] &&
	  pdbRef[i].selected()) {
	Vector dv;

	dv=crd[i]-s.crd[i];

	if (alt[i]>=0) {
	  double d1=dv*dv;
	  dv=crd[i]-s.crd[alt[i]];
	  double d2=dv*dv;
	  
	  if (d1<d2) {
	    d+=d1;
	  } else {
	    d+=d2;
	  }
	} else {
	  d+=dv*dv;
	}
	nd++;
      }
    }
  }

  return sqrt(d/(double)nd);
}

double Structure::rmsd(Structure &s) {
  if (clusterMode == RMSDMODE) {
    return rmsdcoor(s,lsqfitFlag);
  } else if (clusterMode == PHIMODE) {
    return rmsdphi(s);
  } else if (clusterMode == PSIMODE) {
    return rmsdpsi(s);
  } else if (clusterMode == PHIPSIMODE) {
    return 0.5*(rmsdphi(s)+rmsdpsi(s));
  } else if (clusterMode == MIXMODE) {
    return mixfactor*0.05*(rmsdphi(s)+rmsdpsi(s))+
      (1.0-mixfactor)*(rmsdcoor(s,lsqfitFlag));
  } else {
    fprintf(stderr,"unknown cluster mode\n");
    return 0.0;
  }
}

double Structure::dihedral(int i, int j, int k, int l) {
  if (i<0 || j<0 || k<0 || l<0) 
    return 0.0;

  Vector d12=crd[i]-crd[j];
  Vector d23=crd[j]-crd[k];
  Vector d43=crd[l]-crd[k];
 
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

void Structure::calcDihedrals() {
  int *cinx=new int[MAXRES];
  int *ninx=new int[MAXRES];
  int *cainx=new int[MAXRES];

  resLookup=new int[MAXRES];
  nres=0;

  int lastres=-999;

  int i;
  for (i=0; i<natoms; i++) {
    if (pdbRef[i].residueNumber() != lastres) {
      cinx[nres]=ninx[nres]=cainx[nres]=-1;
      resLookup[nres++]=pdbRef[i].residueNumber();
      lastres=pdbRef[i].residueNumber();
      //      fprintf(stderr,"nres: %d, lookup %d -> %d\n",
      //	      nres,nres-1,resLookup[nres-1]);
    }
    if (!strcmp(pdbRef[i].atomName(),"C")) {
      cinx[nres-1]=i;
    } else if (!strcmp(pdbRef[i].atomName(),"N")) {
      ninx[nres-1]=i;
    } else if (!strcmp(pdbRef[i].atomName(),"CA")) {
      cainx[nres-1]=i;
    }
  }

  phi=new double[nres];
  psi=new double[nres];

  for (i=0; i<nres; i++) {
    phi[i]=(i>0)?dihedral(cinx[i-1],ninx[i],cainx[i],cinx[i]):0.0;
    psi[i]=(i<nres-1)?dihedral(ninx[i],cainx[i],cinx[i],ninx[i+1]):0.0;
  }

  delete cinx;
  delete ninx;
  delete cainx;
}

// *** Cluster ****************************************

class Cluster {
 private:
  Structure** sList;
  int nelements;
  int maxelements;

  Structure centr;
  Structure lastcentr;

  double err;

  int havelast;
  
 public:
  Cluster(int maxel);
  ~Cluster();
  
  Structure*& operator[](int inx) { return sList[inx]; }
  int elements() const { return nelements; }

  Structure& centroid()     { return centr; }  
  Structure& lastCentroid() { return lastcentr; } 

  double error() { return err; }

  void calcError();

  void addElement(Structure& s, int lsqfit);
  void clearList();

  void merge(Cluster& c, int lsqfit);

  void copyCentroid();

  int haveLast() { return havelast; }
};
 
Cluster::Cluster(int maxel) : 
  maxelements(maxel), nelements(0), err(-1.0), havelast(0) {
  sList=new Structure*[maxel];
}

Cluster::~Cluster() {
  delete sList;
}

void Cluster::calcError() {
  if (haveLast()) 
    err=centroid().rmsd(lastCentroid());
  else 
    err=-1.0;
}

void Cluster::clearList() {
  nelements=0;
}

void Cluster::merge(Cluster& c, int lsqfit) {
  for (int i=0; i<c.elements(); i++)
    addElement(*c[i],lsqfit);
}

void Cluster::addElement(Structure& s, int lsqfit) {
  if (nelements<maxelements) {
    sList[nelements]=&s;
    if (nelements==0) {
      centr=s;
    } else {
      centr.add(s,nelements,lsqfit);
    }
    nelements++;
  } else {
    fprintf(stderr,"maximum number of elements %d reached\n",maxelements);
  }
}

void Cluster::copyCentroid() {
  lastcentr=centr;
  havelast=1;
}

// *** main ********************************

int main(int argc, char **argv) {
  InputModeEnum inpMode=PDBMODE;
  SelEnum selMode=ALL;
  
  char inputFile[1024];
  *inputFile=0;
  
  double radius=2.5;
  double maxerr=0.5;
  int iterateFlag=0;
  int niter=0;
  
  int centroid=0;
  int cdist=0;
  
  int i,j;

  for (i=0; i<MAXRES; i++) 
    dofit[i]=dormsd[i]=1;

  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i],"-maxerr") && i<argc-1) {
      maxerr=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-mode") && i<argc-1) {
      if (!strcmp(argv[++i],"rmsd"))
	clusterMode=RMSDMODE;
      else if (!strcmp(argv[i],"phipsi")) 
	clusterMode=PHIPSIMODE;
      else if (!strcmp(argv[i],"phi")) 
	clusterMode=PHIMODE;
      else if (!strcmp(argv[i],"psi"))
	clusterMode=PSIMODE;
      else if (!strcmp(argv[i],"mix")) {
	clusterMode=MIXMODE;
	mixfactor=atof(argv[++i]);
      } else
	fprintf(stderr,"unknown cluster mode\n");
    } else if (!strcmp(argv[i],"-radius") && i<argc-1) {
      radius=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-lsqfit")) {
      lsqfitFlag=1;
    } else if (!strcmp(argv[i],"-pdb")) {
      inpMode=PDBMODE;
    } else if (!strcmp(argv[i],"-sicho")) {
      inpMode=SICHOMODE;
    } else if (!strcmp(argv[i],"-iterate")) {
      iterateFlag=1;
    } else if (!strcmp(argv[i],"-centroid")) {
      centroid=1;
    } else if (!strcmp(argv[i],"-cdist")) {
      cdist=1;
    } else if (!strcmp(argv[i],"-ca")) {
      selMode=CA;
    } else if (!strcmp(argv[i],"-cb")) {
      selMode=CB;
    } else if (!strcmp(argv[i],"-cab")) {
      selMode=CAB;
    } else if (!strcmp(argv[i],"-heavy")) {
      selMode=HEAVY;
    } else if (!strcmp(argv[i],"-all")) {
      selMode=ALL;
    } else if (!strcmp(argv[i],"-l")) {
      for (j=0; j<MAXRES; dormsd[j++]=0);
      Field lf(argv[++i],"=");
      for (int ilf=0; ilf<lf.number(); ilf++) {
	Field rf(lf[ilf],":");
	for (j=atoi(rf[0]); j<=atoi(rf[1]); j++) 
	  dormsd[j]=1;
      }
    } else if (!strcmp(argv[i],"-fit")) {
      for (j=0; j<MAXRES; dofit[j++]=0);
      Field lf(argv[++i],"=");
      for (int ilf=0; ilf<lf.number(); ilf++) {
	Field rf(lf[ilf],":");
	for (j=atoi(rf[0]); j<=atoi(rf[1]); j++) 
	  dofit[j]=1;
      }
      lsqfitFlag=1;
    } else if (!strcmp(argv[i],"-fitxl")) {
      for (j=0; j<MAXRES; j++)
	dofit[j]=(dormsd[j])?0:1;
      lsqfitFlag=1;
    } else if (!strcmp(argv[i],"-help")) {
      usage();
    } else if (*argv[i]!='-') {
      strcpy(inputFile,argv[i]);
      break;
    } else {
      fprintf (stderr,"unknown command line option %s\n",argv[i]);
      exit(1);
    }
  }

  if (inpMode == SICHOMODE) 
    clusterMode=RMSDMODE;

  FILE *finp;
  finp=(*inputFile==0) ? stdin : fopen(inputFile,"r");
  if (finp==0) {
    fprintf(stderr,"cannot read input file %s\n",inputFile);
    exit(1);
  }

  Structure **s=new Structure*[MAX_S];
  char **fname=new char*[MAX_S];
  int nstruct=0;

  while (!feof(finp)) {
    if (nstruct>=MAX_S) {
      fprintf(stderr,"Maximum number of input conformations reached (limit is: %d)\n",MAX_S);
      exit(1);
    }
    fname[nstruct]=new char[256];
    if (fscanf(finp,"%s",fname[nstruct])>0) {
      Structure *ts=new Structure(nstruct);
      if (verbose) fprintf(stderr,"reading %d: %s\n",nstruct,fname[nstruct]);
      ts->read(fname[nstruct],inpMode,selMode,pdbRef);
      if (clusterMode!=RMSDMODE) 
	ts->calcDihedrals();
      s[nstruct++]=ts;
    }
  }
  fclose(finp);

  Cluster **cList=new Cluster*[MAX_S];
  int ncluster=0;

  double err=-1.0;
  do {
    for (int is=0; is<nstruct; is++) {
      //            if (verbose) fprintf(stderr,"working on structure %d\n",is);
      double minval=999999.0;
      int mincluster=-1;
      for (int ic=0; ic<ncluster; ic++) {
	double rmsd=(cList[ic]->haveLast())?
	  cList[ic]->lastCentroid().rmsd(*s[is]):
	  cList[ic]->centroid().rmsd(*s[is]);
	//		fprintf(stderr,"   cluster %d has rmsd %lf\n",ic,rmsd);
	if (rmsd<radius && rmsd<minval) {
	  minval=rmsd;
	  mincluster=ic;
	}
      } 
      //            fprintf(stderr," minval: %lf, mincluster: %d\n",minval,mincluster);
      if (mincluster>=0) {
	cList[mincluster]->addElement(*s[is],lsqfitFlag);
	//	if (verbose) fprintf(stderr," adding to cluster %d\n",mincluster);
      } else {
	if (verbose) fprintf(stderr," creating new cluster %d\n",ncluster);
	cList[ncluster]=new Cluster(nstruct);
	cList[ncluster++]->addElement(*s[is],lsqfitFlag);
        if (ncluster==nstruct/2) 
	  fprintf(stderr,"WARNING: reached %d clusters\n",nstruct/2);
      }
    }

    if (iterateFlag) {
      for (i=0; i<ncluster; i++) {
	if (cList[i]->elements() == 0) {
	  //	  fprintf(stderr,"empty cluster %d nuked\n",i);
	  delete cList[i];
	  cList[i]=cList[--ncluster];
	  i--;
	} 
      }

      if (err<0.0) {
	for (i=0; i<ncluster; i++) {
	  int *minc=new int[cList[i]->elements()];
	  int ie=0;
	  do {
	    Structure *s=(*cList[i])[ie];
	    minc[ie]=-1;
	    double minrmsd=99999.0;
	    for (int ic=0; ic<ncluster; ic++) {
	      if (ic!=i) {
		double rmsd=cList[ic]->centroid().rmsd(*s);
		if (rmsd<radius && rmsd<minrmsd) {
		  minrmsd=rmsd;
		  minc[ie]=ic;
		}
	      }
	    }
            if (verbose)
	    	    fprintf(stderr,"cluster %d el %d minc %d minrmsd %lf\n",
	    		    i,ie,minc[ie],minrmsd);
	    ie++;
	  } while(minc[ie-1]>=0 && ie<cList[i]->elements());
	  
	  if (minc[ie-1]>=0) {
	    if (verbose) fprintf(stderr,"elimnating cluster %d\n",i);
	    for (ie=0; ie<cList[i]->elements(); ie++) {
	      Structure *s=(*cList[i])[ie];
	      cList[minc[ie]]->addElement(*s,lsqfitFlag);
	    }
	    delete cList[i];
	    cList[i]=cList[--ncluster];
	    i--;
	  }
	  delete minc;
	}
      }

      /*
      for (i=0; i<ncluster; i++) {
	fprintf(stderr,"cluster %d ",i);
	for (j=0; j<cList[i]->elements(); j++) 
	  fprintf(stderr,"%d ",(*cList[i])[j]->index());
	fprintf(stderr,"\n");
      }
      */

      err=0.0;
      for (i=0; i<ncluster; i++) {
	if (cList[i]->haveLast()) {
	  cList[i]->calcError();
	  //	  fprintf(stderr,"error %d: %lf\n",i,cList[i]->error());
	  err+=cList[i]->error();
	} else {
	  err+=9999.0;
	}
      }
      err/=(double)ncluster;

      if (++niter>10) 
	iterateFlag=0;

      if (verbose) fprintf(stderr,"iteration %d: error = %lf\n",niter,err);

      for (i=0; i<ncluster; i++) {
	cList[i]->copyCentroid();
	
	if (err>maxerr && iterateFlag) 
	  cList[i]->clearList();
      }
    }
  } while (iterateFlag && err>maxerr);

  /*  
  for (i=0; i<ncluster; i++) {
    fprintf(stderr,"cluster %d\n",i);
    for (int ie=0; ie<cList[i]->elements(); ie++) {
      Structure *s=(*cList[i])[ie];
      fprintf(stderr," %2d ",s->index());
      for (int ic=0; ic<ncluster; ic++) {
	double rmsd=cList[ic]->centroid().rmsd(*s);
	fprintf(stderr," %lf",rmsd);
      }
      fprintf(stderr,"\n");
    }
  }
  */
  
  for (i=0; i<ncluster; i++) {
    fprintf(stdout,"#Cluster %d\n",i+1);
    for (int ie=0; ie<cList[i]->elements(); ie++) {
      Structure *ts=(*cList[i])[ie];
      int inx=ts->index();
      if (cdist) {
	fprintf(stdout,"%d %s %lf\n",inx+1,fname[inx],cList[i]->centroid().rmsd(*ts));
      } else {
	fprintf(stdout,"%d %s\n",inx+1,fname[inx]);
      }
    }

    if (centroid && (clusterMode==RMSDMODE || clusterMode==MIXMODE)) {
      fprintf(stdout,"#Centroid %d\n",i+1);
      for (j=0; j<cList[i]->centroid().atoms(); j++) {
	if (inpMode==PDBMODE) {
	  PDBEntry pdb(j+1,pdbRef[j].atomName(),pdbRef[j].residueName(),
		       pdbRef[j].residueNumber(),cList[i]->centroid().coordinates(j));
	  pdb.write(stdout);
	} else {
	  fprintf(stdout,"%lf %lf %lf\n",
		  cList[i]->centroid().coordinates(j).x(),
		  cList[i]->centroid().coordinates(j).y(),
		  cList[i]->centroid().coordinates(j).z());
	}
      }
    }
    fprintf(stdout,"#End\n");
  }

  for (i=0; i<ncluster; i++) 
    delete cList[i];
  delete cList;

  delete fname;
  delete s;
}
  
  

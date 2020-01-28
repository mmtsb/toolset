// http://mmtsb.scripps.edu/doc/jclust.html
// 2000, Michael Feig, Brooks group, TSRI
// 2000, John Karanicolas, Brooks group, TSRI
//
// implementing algorithms from the following papers:
//
// T. Kurita: An Efficient Agglomerative Clustering Algorithm
//   Using a Heap, Pattern Recognition (1991), vol. 24, pp. 204-209
// S. Xu, M. V. Kamath, D. W. Capson: Selection of Partitions 
//   from a Hierarchy, Pattern Recognition Letters (1993), vol. 14,
//   pp. 7-15
//
// original code by John Karanicolas
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
  fprintf(stderr,"usage:   jclust [options] [inputFile]\n");
  fprintf(stderr,"options: [-mode rmsd|phipsi|psi|phi|contact]\n");
  fprintf(stderr,"         [-maxdist value]\n");
  fprintf(stderr,"         [-centroid] [-ref file]\n");
  fprintf(stderr,"         [-pdb | -sicho] [-max value]\n");
  fprintf(stderr,"         [-ca | -cb | -cab | -heavy | -all]\n");
  fprintf(stderr,"         [-l min:max[=min:max=...]]\n");
  fprintf(stderr,"         [-lsqfit] [-fitxl] [-fit min:max[=min:max=...]]\n");
  exit(1);
}

enum InputModeEnum { PDBMODE, SICHOMODE };
enum ClusterModeEnum { RMSDMODE, CONTACTMODE, PHIPSIMODE, PHIMODE, PSIMODE };

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
const int MAXRES=4096;

const int MAX_S=32000;

int dofit[MAXRES];
int dormsd[MAXRES];

PDBEntry pdbRef[MAXATOMS];

// *** Structure ***********************************

class Structure {
 private:
  Vector *crd;
  double *phi;
  double *psi;

  int nres;
  int natoms;
  int inx;

  int *resLookup;
  
 public:
  Structure(int i=0);
  Structure(const Structure& conf);
  ~Structure();
  
  Structure& operator=(const Structure& conf);
  Vector& operator[](int i) { return crd[i]; }  
  
  Vector& coordinates(int inx) { return crd[inx]; }
  double& phiDihedral(int inx) { return phi[inx]; }
  double& psiDihedral(int inx) { return psi[inx]; }

  int residues() const { return nres; }
  int atoms() const { return natoms; }
  int index() const { return inx; }
  
  void read(char *fname, InputModeEnum=PDBMODE, 
	    SelEnum sel=ALL, PDBEntry *pdb=0);

  void calcDihedrals();

  double rmsd(Structure& s, int lsqFlag);

  double rmsdcoor(Structure& s, int lsqFlag);
  double rmsdphi(Structure& s);
  double rmsdpsi(Structure& s);

  double rmsd(Structure& s, ClusterModeEnum mode, int lsqFlag=1);

  void lsqFit(Structure& s);

  double dihedral(int i, int j, int k, int l);
};

Structure::Structure(int i) : 
  natoms(0), inx(i), crd(0), phi(0), psi(0), resLookup(0), nres(0) {}

Structure::~Structure() {
  if (crd!=0) 
    delete[] crd;

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

  int i;
  crd=new Vector[natoms];
  for (i=0; i<natoms; i++)
    crd[i]=conf.crd[i];

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
	if (pdb.selected()) {
	  if (pdbs[n].atomIndex()<=0)
	    pdbs[n]=pdb;
	  v[n++]=pdb.coordinates();
	}
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
  } else if (natoms!=n) {
    fprintf(stderr,"number of atoms do not match (file: %d, init: %d)\n",n,natoms);
    exit(1);
  }
  
  for (int i=0; i<natoms; i++)
    crd[i]=v[i];

  delete[] v;
}

void Structure::lsqFit(Structure& s) {
  Vector tc1(0,0,0);
  Vector tc2(0,0,0);
    
  int i;
  int nlsqfit=0;
  for (i=0; i<natoms; i++) {
    if (dofit[pdbRef[i].residueNumber()]) {
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
    if (dofit[pdbRef[i].residueNumber()]) {
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
  for (i=1; i<nres-1; i++) {
    if (dormsd[resLookup[i]]) {
      double dv=(phi[i]-s.phi[i]);
      if (dv>180.0) 
	dv-=360.0;
      else if (dv<-180.0) 
	dv+=360.0;
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

  for (i=1; i<nres-1; i++) {
    if (dormsd[resLookup[i]]) {
      double dv=psi[i]-s.psi[i];
      if (dv>180.0) 
	dv-=360.0;
      else if (dv<-180.0) 
	dv+=360.0;
      d+=dv*dv;
      nd++;
    }
  }
  return sqrt(d/(double)nd);
}

double Structure::rmsdcoor(Structure& s, int lsqFlag) {
  register int i;
  register double d=0.0;

  if (lsqFlag) {
    Structure ns(s);
    ns.lsqFit(*this);
    for (i=0; i<natoms; i++) {
      if (dormsd[pdbRef[i].residueNumber()]) {
	Vector dv;
	dv=crd[i]-ns.crd[i];
	d+=dv*dv;
      }
    }
  } else {
    for (i=0; i<natoms; i++) {
      if (dormsd[pdbRef[i].residueNumber()]) {
	Vector dv;
	dv=crd[i]-s.crd[i];
	d+=dv*dv;
      }
    }
  }

  return d;
}

double Structure::rmsd(Structure &s, ClusterModeEnum clusterMode, int lsqFlag) {
  if (clusterMode == RMSDMODE) {
    return rmsdcoor(s,lsqFlag);
  } else if (clusterMode == PHIMODE) {
    return rmsdphi(s);
  } else if (clusterMode == PSIMODE) {
    return rmsdpsi(s);
  } else if (clusterMode == PHIPSIMODE) {
    return 0.5*(rmsdphi(s)+rmsdpsi(s));
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

// *** ContactMap **********************************

class ContactMap {
 private:
  double maxdist;

  int structInx;
  
  int natoms;

  double *dist;
  int ndist;
  int naccu;

  int readFile(char *fname, InputModeEnum mode, Vector *crd);

 public:
  ContactMap();
  ContactMap(int inx, char *fname, InputModeEnum mode, double mxd=12.0);
  ~ContactMap();

  ContactMap& operator=(const ContactMap&);

  double compare(ContactMap& cmap);
  void accumulate(ContactMap& cmap);

  int elements() const { return naccu; }
  int index() const { return structInx; }

  void writeRGB(FILE *fptr, ContactMap *refMap=0);
};

void ContactMap::writeRGB(FILE *fptr, ContactMap *refMap) {
  int *r=new int[natoms*natoms];
  int *g=new int[natoms*natoms];
  int *b=new int[natoms*natoms];

  int i,j;

  for (i=0; i<natoms*natoms; i++) {
    r[i]=g[i]=b[i]=0;
  }

  int rval[256];
  int gval[256];
  int bval[256];

  for (i=0; i<256; i++) {
    double tv=(double)i/180.0*M_PI;
    rval[i]=int((cos(tv-120.0/180.0*M_PI)*127.0))+128;
    gval[i]=int((cos(tv)*127.0))+128;
    bval[i]=int((cos(tv+120.0/180.0*M_PI)*127.0))+128;
  }

  int nn=0;
  for (i=0; i<natoms-4; i++) {
    for (j=i+4; j<natoms; j++) {
      int inx=i*natoms+j;
      int val=int((dist[nn++]-3)*28.0);
      if (val>255) 
	val=255;
      if (val<0) 
	val=0;
      r[inx]=rval[val];
      g[inx]=gval[val];
      b[inx]=bval[val];
    }
  }

  if (refMap!=0) {
    nn=0;
    for (i=0; i<natoms-4; i++) {
      for (j=i+4; j<natoms; j++) {
	int inx=j*natoms+i;
	int val=int((refMap->dist[nn++]-3)*28.0);
	if (val>255) 
	  val=255;
	if (val<0) 
	  val=0;
	r[inx]=rval[val];
	g[inx]=gval[val];
	b[inx]=bval[val];
      }
    }
  }

  /*
  for (j=0; j<natoms; j++) {
    int inx=(natoms-1)*natoms+j;
    int val=int((double)j/(double)natoms*256.0);
    r[inx]=rval[val];
    g[inx]=gval[val];
    b[inx]=bval[val];
  }
  */
  
  fprintf(fptr,"P3\n%3d %3d\n256",natoms,natoms);
  for (i=0; i<natoms*natoms; i++) {
    if (i%5==0) 
      fprintf(fptr,"\n");
    fprintf(fptr,"%3d %3d %3d ",r[i],g[i],b[i]);
  }  
  fprintf(fptr,"\n");
  
  delete r;
  delete g;
  delete b;
}

ContactMap& ContactMap::operator=(const ContactMap& from) {
  maxdist=from.maxdist;
  structInx=from.structInx;
  ndist=from.ndist;
  naccu=from.naccu;
  natoms=from.natoms;

  if (dist!=0) 
    delete dist;

  if (ndist>0) {
    dist=new double[ndist];
    for (int i=0; i<ndist; i++) 
      dist[i]=from.dist[i];
  }

  return *this;
}

ContactMap::ContactMap() :
  structInx(-1), ndist(0), naccu(0), dist(0), natoms(0), maxdist(12.0) {}

ContactMap::ContactMap(int inx, char *fname, InputModeEnum mode, double mxd) : 
  structInx(inx), ndist(0), naccu(1), maxdist(mxd) {
  Vector *crd=new Vector[5000];
  natoms=readFile(fname,mode,crd);
  dist=new double[natoms*(natoms-1)/2];
  for (int i=0; i<natoms-4; i++) {
    for (int j=i+4; j<natoms; j++) {
      Vector v=crd[i]-crd[j];
      dist[ndist++]=v.norm();
    }
  }
  delete[] crd;
}

ContactMap::~ContactMap() {
  if (dist!=0) 
    delete dist;
}

double ContactMap::compare(ContactMap& cmap) {
  double avg=0.0;
  //  int navg=0;
  for (int i=0; i<ndist; i++) {
    double v=cmap.dist[i];
    if (v>maxdist) 
      v=maxdist;

    double w=dist[i];
    if (w>maxdist) 
      w=maxdist;

    double d=v-w;
    avg+=d*d;
  }
  avg/=(double)ndist;
  return sqrt(avg);
}

void ContactMap::accumulate(ContactMap& cmap) {
  for (int i=0; i<ndist; i++) {
    dist[i]=(dist[i]*(double)naccu+cmap.dist[i]*(double)cmap.naccu)/
      (double)(naccu+cmap.naccu);
  }
  naccu++;
}

int ContactMap::readFile(char *fname, InputModeEnum mode, Vector *crd) {
  static int lastN=0;

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

  int num;
  int n=0;
  while (!feof(fptr)) {
    if (mode == PDBMODE) {
      if (pdb.read(fptr,CB)>0) {
	if (pdb.selected()) {
	  crd[n]=pdb.coordinates();
          n++;
	}
      }
    } else if (mode == SICHOMODE) {
      if (n==0) {
	fscanf(fptr,"%d",&num);
	n++;
      } else {
	int ix,iy,iz;
	if (fscanf(fptr,"%5d%5d%5d",&ix,&iy,&iz)>0) {
          if (n>1 && n<num) 
	    crd[n-2]=Vector((double)ix,(double)iy,(double)iz);
	  n++;
	} 
      }
    }
  }
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

  return n;
}

// *** Link *****************************

class Link {
 private:
  void* el;
  Link* nxt;

 public:
  Link(void* s, Link* n=NULL);
  ~Link();
  
  void*& element() { return el; }
  Link*& next()  { return nxt; }
};

Link::Link(void *e, Link *n) :
  el(e), nxt(n) {}

Link::~Link() {}

// *** Cluster ***********************************

class Cluster {
 protected:
  Link* sList;
  Link* last;
  int nelements;
  double sse;

  void clearList();

 public:
  Cluster();
  ~Cluster();

  Link*& elementList()  { return sList; }
  int elements() const  { return nelements; }

  void addElement (void *s);
  virtual void setCentroid (void *s)=0;

  virtual double merge(Cluster* c, int lsqFlag=1)=0;
  virtual double compare(Cluster* c, int lsqFlag=1)=0;
};

Cluster::Cluster() : 
  sList(NULL), last(NULL), nelements(0), sse(0.0) {
}

Cluster::~Cluster() {
  clearList();
}
    
void Cluster::clearList() {
  Link *l=sList;
  while (l!=NULL) {
    Link *t=l;
    l=l->next();
    delete t;
  }
  sList=NULL;
}

void Cluster::addElement(void *s) {
  Link *newLink=new Link(s,sList);
  sList=newLink;
  if (last==NULL) 
    last=sList;
  nelements++;
}

// *** RMSDCluster ***********************************

class RMSDCluster : public Cluster {
 private:
  Structure centr;
  ClusterModeEnum clusterMode;

 public:
  RMSDCluster(ClusterModeEnum cMode);
  RMSDCluster& operator=(const RMSDCluster & c);
  Structure& centroid() { return centr; }
  double merge(Cluster* c, int lsqFlag=1);
  void setCentroid (void *s);
  double compare(Cluster* c, int lsqFlag=1);
};

RMSDCluster::RMSDCluster(ClusterModeEnum cMode) : 
  Cluster(), centr(-1), clusterMode(cMode) {}

void RMSDCluster::setCentroid(void *s) {
  centr=*(Structure *)s;
}

double RMSDCluster::compare(Cluster *c, int lsqFlag) {
  RMSDCluster *rc=(RMSDCluster *)c;
  return centr.rmsd(rc->centr,clusterMode,lsqFlag);
}

double RMSDCluster::merge(Cluster* c, int lsqFlag) {
  RMSDCluster *rc=(RMSDCluster *)c;
  last->next()=rc->sList;
  last=rc->last;
  rc->sList=NULL;

  double csse=-rc->sse-sse;

  if (lsqFlag && clusterMode==RMSDMODE) 
    rc->centr.lsqFit(centr);

  int comb=rc->nelements+nelements;
  
  if (clusterMode==RMSDMODE) {
    for (int i=0; i<centr.atoms(); i++) 
      centr[i]=(centr[i]*nelements+(rc->centr)[i]*rc->nelements)/(double)comb;
  } else {
    for (int i=0; i<centr.residues(); i++) {
      double tphi=centr.phiDihedral(i);
      double rphi=rc->centr.phiDihedral(i);
      if (rphi>tphi+180.0)
        rphi-=360.0;
      else if (rphi<tphi-180.0) 
        rphi+=360.0;
      centr.phiDihedral(i)=(tphi*nelements+rphi*rc->nelements)/(double)comb;
      if (centr.phiDihedral(i)<-180.0) 
	centr.phiDihedral(i)+=360.0;
      else if (centr.phiDihedral(i)>180.0) 
	centr.phiDihedral(i)-=360.0;

      double tpsi=centr.psiDihedral(i);
      double rpsi=rc->centr.psiDihedral(i);
      if (rpsi>tpsi+180.0)
        rpsi-=360.0;
      else if (rpsi<tpsi-180.0) 
        rpsi+=360.0;
      centr.psiDihedral(i)=(tpsi*nelements+rpsi*rc->nelements)/(double)comb;
      if (centr.psiDihedral(i)<-180.0) 
	centr.psiDihedral(i)+=360.0;
      else if (centr.psiDihedral(i)>180.0) 
	centr.psiDihedral(i)-=360.0;
    }
  }

  nelements=comb;
  
  sse=0.0;
  for (Link *l=sList; l!=NULL; l=l->next()) 
    sse+=centr.rmsd(*(Structure *)(l->element()),clusterMode,lsqFlag);

  rc->nelements=0;

  return csse+sse;
}

RMSDCluster& RMSDCluster::operator=(const RMSDCluster& c) {
  centr=c.centr;
  nelements=0;
  sse=c.sse;
  clusterMode=c.clusterMode;
  last=NULL;
  clearList();
  for (Link *l=c.sList; l!=NULL; l=l->next()) {
    addElement(l->element());
  }
  return (*this);
}

// *** ContactCluster ***********************************

class ContactCluster : public Cluster {
 private:
  ContactMap centr;

 public:
  ContactCluster();
  ContactCluster& operator=(const ContactCluster & c);
  ContactMap& centroid() { return centr; }
  double merge(Cluster* c, int);
  void setCentroid(void *s);
  double compare(Cluster* c, int);
};

ContactCluster::ContactCluster() : 
  Cluster(), centr() {}

void ContactCluster::setCentroid(void *s) {
  centr=*(ContactMap *)s;
}

double ContactCluster::compare(Cluster *c, int) {
  ContactCluster *cc=(ContactCluster *)c;
  return centr.compare(cc->centr);
}

double ContactCluster::merge(Cluster* c, int) {
  ContactCluster *cc=(ContactCluster *)c;

  last->next()=cc->sList;
  last=cc->last;
  cc->sList=NULL;

  double csse=-cc->sse-sse;

  centr.accumulate(cc->centr);
  nelements+=cc->nelements;
  
  sse=0.0;
  for (Link *l=sList; l!=NULL; l=l->next()) 
    sse+=centr.compare(*(ContactMap *)(l->element()));

  cc->nelements=0;

  return csse+sse;
}


ContactCluster& ContactCluster::operator=(const ContactCluster& c) {
  centr=c.centr;
  nelements=0;
  sse=c.sse;
  last=NULL;
  clearList();
  for (Link *l=c.sList; l!=NULL; l=l->next()) {
    addElement(l->element());
  }
  return (*this);
}

// *** global functions ******

void shiftdown(double *d, int *p, int *w, int s, int e) {
  int i,j;
  double x;
  int ip,ipp;
  
  i=s; j=2*i+1;
  x=d[i];ip=p[i];
  while(j<=e) {
    if (j<e) {
      if (d[j]>d[j+1]) 
	j++;
    }
    if (x<=d[j]) break;
    d[i]=d[j]; p[i]=ipp=p[j];
    w[ipp]=i;
    i=j; j=2*i+1;
  }
  d[i]=x; w[ip]=i; p[i]=ip;
}

void shiftup(double *d, int *p, int *w, int s) {
  int i,j;
  double x;
  int ip,ipp;

  i=s; j=((i+1)/2)-1; 
  x=d[i]; ip=p[i];
  while(j>=0) {
    if (d[j]<=x) break;
    d[i]=d[j]; p[i]=ipp=p[j];
    w[ipp]=i;
    i=j; j=((i+1)/2)-1;
  }
  d[i]=x; p[i]=ip; w[ip]=i;
}

// *** main ********************************

int main(int argc, char **argv) {
  int maxnofcluster=-1;

  InputModeEnum inpMode=PDBMODE;
  ClusterModeEnum clusterMode=RMSDMODE; 
  SelEnum selMode=ALL;
  
  int lsqfitFlag=0;

  char inputFile[1024];
  *inputFile=0;

  int centroid=0;

  char refFile[1024];
  *refFile=0;

  double contMaxDist=12.0;

  int i,j;

  for (i=0; i<MAXRES; i++) 
    dofit[i]=dormsd[i]=1;

  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i],"-max") && i<argc-1) {
      maxnofcluster=atoi(argv[++i]);
    } else if (!strcmp(argv[i],"-mode") && i<argc-1) {
      if (!strcmp(argv[++i],"rmsd"))
	clusterMode=RMSDMODE;
      else if (!strcmp(argv[i],"contact")) 
	clusterMode=CONTACTMODE;
      else if (!strcmp(argv[i],"phi")) 
	clusterMode=PHIMODE;
      else if (!strcmp(argv[i],"psi"))
	clusterMode=PSIMODE;
      else if (!strcmp(argv[i],"phipsi"))
	clusterMode=PHIPSIMODE;
      else 
	fprintf(stderr,"unknown cluster mode\n");
    } else if (!strcmp(argv[i],"-maxdist")) {
      contMaxDist=atoi(argv[++i]);
    } else if (!strcmp(argv[i],"-lsqfit")) {
      lsqfitFlag=1;
    } else if (!strcmp(argv[i],"-pdb")) {
      inpMode=PDBMODE;
    } else if (!strcmp(argv[i],"-sicho")) {
      inpMode=SICHOMODE;
    } else if (!strcmp(argv[i],"-centroid")) {
      centroid=1;
    } else if (!strcmp(argv[i],"-ref")) {
      strcpy(refFile,argv[++i]);
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

  FILE *finp;
  finp=(*inputFile==0) ? stdin : fopen(inputFile,"r");
  if (finp==0) {
    fprintf(stderr,"cannot read input file %s\n",inputFile);
    exit(1);
  }

  void *s[MAX_S];
  char *fname[MAX_S];

  int nstruct=0;

  while (!feof(finp)) {
    if (nstruct>=MAX_S) {
      fprintf(stderr,"Maximum number of input conformations reached (limit is: %d)\n",MAX_S);
      exit(1);
    }
    fname[nstruct]=new char[256];
    if (fscanf(finp,"%s",fname[nstruct])>0) {
      if (clusterMode!=CONTACTMODE) {
	//	fprintf(stderr,"read file %s\n",fname[nstruct]);
	Structure *ts=new Structure(nstruct);
	ts->read(fname[nstruct],inpMode,selMode,pdbRef);
	if (clusterMode!=RMSDMODE) 
	  ts->calcDihedrals();
	s[nstruct++]=(void *)ts;
      } else {
	ContactMap *ts=new ContactMap(nstruct,fname[nstruct],
				      inpMode,contMaxDist);
	s[nstruct++]=(void *)ts;
      }
    }
  }
  fclose(finp);

  fprintf(stderr,"running jclust (mode %d) for %d structures\n",
	  clusterMode,nstruct);

  int maxclust=(maxnofcluster<0)?nstruct/2:maxnofcluster;

  int nrClust;  
  int* valid=new int[nstruct];

  Cluster** rClust=new Cluster*[nstruct];
  Cluster** clust=new Cluster*[nstruct];

  for (i=0; i<nstruct; i++) {
    if (clusterMode!=CONTACTMODE) {
      rClust[i]=new RMSDCluster(clusterMode);
      clust[i]=new RMSDCluster(clusterMode);
    } else {
      rClust[i]=new ContactCluster;
      clust[i]=new ContactCluster;
    }
    clust[i]->setCentroid(s[i]);
    clust[i]->addElement(s[i]);
    valid[i]=1;
  }

  int m=(nstruct*(nstruct-1)/2);

  double* cDist= new double[m];
  int* cPair= new int[m];
  int* cWhere= new int[m];
  
  int *iLookup= new int[nstruct*(nstruct-1)/2];
  int *jLookup= new int[nstruct*(nstruct-1)/2];
  int *nLookup= new int[nstruct*nstruct];

  int n=0;
  for (i=0; i<nstruct-1; i++) {
    //    fprintf(stderr,"compare distances for %d\n",i);
    for (j=i+1; j<nstruct; j++) {
      cDist[n]=0.5*clust[i]->compare(clust[j],lsqfitFlag);
      // fprintf(stderr,"compare %d <-> %d : %lf\n",i,j,cDist[n]);
      cPair[n]=n;
      cWhere[n]=n;
      iLookup[n]=i;
      jLookup[n]=j;
      nLookup[i*nstruct+j]=n;
      n++;
    }
  }

  int tn=m/2; 
  while(tn>0) 
    shiftdown(cDist,cPair,cWhere,--tn,m-1);

  double* totalSSE=new double[nstruct];
  double* mbcd=new double[nstruct];

  double maxxif=0.0;

  totalSSE[0]=0.0;
  mbcd[0]=sqrt(cDist[0]);

  for (int k=1; k<nstruct-1; k++) {
    //    fprintf(stderr,"process cluster %d\n",k);
    int isw=iLookup[cPair[0]];
    int jsw=jLookup[cPair[0]];

    totalSSE[k]=totalSSE[k-1]+clust[isw]->merge(clust[jsw],lsqfitFlag);
    valid[jsw]=0;

    for (i=0; i<nstruct; i++) {
      if (i!=isw && valid[i]) {
	n=cWhere[nLookup[(i<isw)?(i*nstruct+isw):(isw*nstruct+i)]];
	
	double newval=clust[isw]->compare(clust[i],lsqfitFlag);
	double weight=(double)(clust[isw]->elements()*clust[i]->elements())/
		  (double)(clust[isw]->elements()+clust[i]->elements());
	newval*=weight;

	if (newval>cDist[n]) {
	  cDist[n]=newval;
	  shiftdown(cDist,cPair,cWhere,n,m-1);
	} else if (newval<cDist[n]) {
	  cDist[n]=newval;
	  shiftup(cDist,cPair,cWhere,n);
	}
      }
    }

    for (i=0; i<nstruct; i++) {
      if (valid[i]) {
	n=cWhere[nLookup[(i<jsw)?(i*nstruct+jsw):(jsw*nstruct+i)]];
	
	int ip;
	double xval=cDist[n];
	cDist[n]=cDist[m-1];
	cPair[n]=ip=cPair[m-1];
	cWhere[ip]=n;
	if (xval<cDist[m-1]) {
	  shiftdown(cDist,cPair,cWhere,n,m-2);
	} else if (xval>cDist[m-1]) {
	  shiftup(cDist,cPair,cWhere,n);
	}
	m--;
      }
    }

    mbcd[k]=sqrt(cDist[0]);

    int nval=0;
    for (int ival=0; ival<nstruct; ival++) 
      if (valid[ival]) nval++;

    double deltaSSE=sqrt(totalSSE[k])-sqrt(totalSSE[k-1]);
    double deltaMBCD=mbcd[k]-mbcd[k-1];
    double xif=(deltaSSE<0.0000001)?-1E99:deltaMBCD/deltaSSE;

    if (deltaSSE<-0.00001 || deltaMBCD<-0.00001) {
      fprintf(stderr,"negative delta: SSE: %lf, MBCD: %lf\n",deltaSSE,deltaMBCD);
      //      exit(1);
    } else {
//      fprintf(stderr,"SSE: %lf\n",deltaSSE);
    }

    if (xif>maxxif && (nstruct-k)<=maxclust) {
      maxxif=xif;
      nrClust=0;
      for (i=0; i<nstruct; i++) {
	if (valid[i]) {
	  if (clusterMode!=CONTACTMODE) {
	    *(RMSDCluster *)rClust[nrClust++]=*(RMSDCluster *)clust[i];
	  } else {
	    *(ContactCluster *)rClust[nrClust++]=*(ContactCluster *)clust[i];
	  }
	}
      }
    }
  }

  for (i=0; i<nrClust; i++) {
    fprintf(stdout,"#Cluster %d\n",i+1);
    for (Link *l=rClust[i]->elementList(); l!=NULL; l=l->next()) {
      int inx;
      if (clusterMode==CONTACTMODE) {
	ContactMap *s=(ContactMap *)l->element();	
	inx=s->index();
      } else {
	Structure *s=(Structure *)l->element();
	inx=s->index();
      }
      fprintf(stdout,"%d %s\n",inx+1,fname[inx]);
    }

    if (centroid) {
      fprintf(stdout,"#Centroid %d\n",i+1);
      if (clusterMode==RMSDMODE) {
	for (j=0; j<((RMSDCluster *)rClust[i])->centroid().atoms(); j++) {
	  if (inpMode==PDBMODE) {
	    PDBEntry pdb(j+1,pdbRef[j].atomName(),pdbRef[j].residueName(),
			 pdbRef[j].residueNumber(),
			 (((RMSDCluster *)rClust[i])->centroid())[j]);
	    pdb.write(stdout);
	  } else {
	    fprintf(stdout,"%lf %lf %lf\n",
		    (((RMSDCluster *)rClust[i])->centroid())[j].x(),
		    (((RMSDCluster *)rClust[i])->centroid())[j].y(),
		    (((RMSDCluster *)rClust[i])->centroid())[j].z());
	  }
	}
      } else if (clusterMode==PHIMODE) {
	for (j=0; j<((RMSDCluster *)rClust[i])->centroid().residues(); j++) {
	  fprintf(stdout,"%lf\n",((RMSDCluster *)rClust[i])->centroid().phiDihedral(j));
	}
      } else if (clusterMode==PSIMODE) {
	for (j=0; j<((RMSDCluster *)rClust[i])->centroid().residues(); j++) {
	  fprintf(stdout,"%lf\n",((RMSDCluster *)rClust[i])->centroid().psiDihedral(j));
	}
      } else if (clusterMode==PHIPSIMODE) {
	for (j=0; j<((RMSDCluster *)rClust[i])->centroid().residues(); j++) {
	  fprintf(stdout,"%lf %lf\n",
		  ((RMSDCluster *)rClust[i])->centroid().phiDihedral(j),
		  ((RMSDCluster *)rClust[i])->centroid().psiDihedral(j));
	}
      } else if (clusterMode==CONTACTMODE) {
	ContactMap *refMap=0;
	if (refFile!=0 && *refFile!=0) {
	  refMap=new ContactMap(0,refFile,inpMode,contMaxDist);
	}

	((ContactCluster *)rClust[i])->centroid().writeRGB(stdout,refMap);
      }
    }

    fprintf(stdout,"#End\n");
  }

  delete valid;
  delete cDist;
  delete cPair;
  delete cWhere;
  delete iLookup;
  delete jLookup;
  delete nLookup;
  delete totalSSE;
  delete mbcd;

  for (i=0; i<nstruct; i++) {
    delete rClust[i];
    delete clust[i];
    if (clusterMode!=CONTACTMODE) {
      delete (Structure *)s[i];
    } else {
      delete (ContactMap *)s[i];
    }
    delete fname[i];
  }
  delete clust;
  delete rClust;
}
  
  

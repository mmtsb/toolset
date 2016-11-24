// http://mmtsb.scripps.edu/doc/rebuild.html
// 2000, Michael Feig, Brooks group, TSRI

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
  fprintf(stderr,"usage:    rebuild [options] [inputfile]\n");
  fprintf(stderr,"options:  [-l min:max[=min:max ...] refpdb]\n");
  fprintf(stderr,"          [-d datadir]\n");
  fprintf(stderr,"          [-backonly] [-fixca]\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"Input file format:\n");
  fprintf(stderr," Index Residue SICHO-x SICHO-y SICHO-z CA-x CA-y CA-z\n");
  fprintf(stderr," ...\n");
  exit(1);
}

#define SCWRLLIB "scwrlbin34.red.lib"
const int MAXRES=5000;
int valid[MAXRES];

AA *res[MAXRES];
int nres=0;

int fixCA;

class ReducedLibRecord {
 public:
  int aatype;
  int phipsi;
  float p;
  float chi1;
  float chi2;
  float chi3;
  float chi4;
};

static int nrotamers[20] = { 81, 9, 9, 3,27, 
                             27, 6, 9, 9,81,
                             27, 6, 2, 3, 3,
 			      9, 6, 3, 0, 0 };

class RotamerData {
 public:
  double p;
  double chi1;
  double chi2;
  double chi3;
  double chi4;
};


class Rotamer {
 private:
  int inum;
  int jnum;
  int *knum;
  int *offset;
  RotamerData *dat;
  int maxdat;
  
  int getIndex(int i, int j, int k);

 public:
  Rotamer(int i, int j, int *k);
  ~Rotamer();

  double chi1(int i, int j, int k);
  double chi2(int i, int j, int k);
  double chi3(int i, int j, int k);
  double chi4(int i, int j, int k);

  void setRotamer(int i, int j, int k, ReducedLibRecord &red);

  int kRotamers(int i);
};

Rotamer::Rotamer(int i, int j, int *k) : inum(i), jnum(j) {
  knum=new int[i];
  offset=new int[i+1];
  offset[0]=0;
  for (int ii=0; ii<i; ii++) {
    knum[ii]=k[ii];
    offset[ii+1]=offset[ii]+j*k[ii];
  }
  maxdat=offset[i];
  dat=new RotamerData[maxdat];
}

Rotamer::~Rotamer() {
  delete[] dat;
  delete knum;
  delete offset;
}

int Rotamer::getIndex(int i, int j, int k) {
  int inx=offset[i]+k*jnum+j;
  if (inx>=maxdat || inx<0) {
    fprintf(stderr,"index %d out of range ( 0:%d )\n",inx,maxdat);
    exit(1);
  }
  return inx;
}

void Rotamer::setRotamer(int i, int j, int k, ReducedLibRecord &red) {
  int inx=getIndex(i,j,k);
  dat[inx].p=(double)red.p;
  dat[inx].chi1=(double)red.chi1;
  dat[inx].chi2=(double)red.chi2;
  dat[inx].chi3=(double)red.chi3;
  dat[inx].chi4=(double)red.chi4;
}

double Rotamer::chi1(int i, int j, int k) {
  return dat[getIndex(i,j,k)].chi1;
}

double Rotamer::chi2(int i, int j, int k) {
  return dat[getIndex(i,j,k)].chi2;
}

double Rotamer::chi3(int i, int j, int k) {
  return dat[getIndex(i,j,k)].chi3;
}

double Rotamer::chi4(int i, int j, int k) {
  return dat[getIndex(i,j,k)].chi4;
}

int Rotamer::kRotamers(int i) {
  return knum[i];
}

int getPhiPsiIndex(double phi, double psi) {
  if (phi<-180.0) 
    phi+=360.0;
  if (phi>=180.0) 
    phi-=360.0;

  if (psi<-180.0) 
    psi+=360.0;
  if (psi>=180.0) 
    psi-=360.0;
  
  int iphi=int((phi/10.0+18));
  int ipsi=int((psi/10.0+18));
  
  return iphi*36+ipsi;
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


int vround(double v) {
  return  (v>0.0)?int(v+0.5):int(v-0.5);
}

int readInputFile(char *fname, AA **res, int *resinx) {
  FILE *fptr;
  fptr= (!strcmp(fname,"-") || *fname==0) ? stdin : fopen(fname,"r");

  if (fptr==0) {
    fprintf(stderr,"cannot open input file %s\n",fname);
    exit(1);
  }
  
  char line[2048];

  int nres=1;

  while (!feof(fptr)) {
    if (fgets(line,2047,fptr)!=0) {
      Field f(line," \n");
      if (f.number()>=5) {
	resinx[atoi(f[0])+100]=nres;
	res[nres]=getAA(f[1]);
	res[nres]->pdbnum=atoi(f[0]);
	res[nres]->sicho=Vector(atof(f[2]),atof(f[3]),atof(f[4]));

	if (f.number()>=8) {
	  res[nres]->modelca=1;
	  res[nres]->ca=Vector(atof(f[5]),atof(f[6]),atof(f[7]));

          // rebuild expects side chain center to include C-alpha
	  //	  int nheavy=res[nres]->heavyAtoms();
	  //	  res[nres]->sicho=
	  //	    (res[nres]->sicho*(double)nheavy+res[nres]->ca)/
	  //	    (double)(nheavy+1);
	} 

	if (++nres>MAXRES-2) {
	  fprintf(stderr,"increase MAXRES\n");
	  exit(1);
	}
      }
    }
  }
  fclose(fptr);
  
  res[0]=getAA("GLY");
  res[nres]=getAA("GLY");
  res[nres+1]=getAA("GLY");

  if (nres<3) {
    fprintf(stderr,"not enough residues\n");
    exit(1);
  }

  Vector d=res[1]->sicho-res[2]->sicho;
  d/=d.norm();
  d.x()+=0.2;
  res[0]->sicho=res[1]->sicho+d*4.5;
  d=res[nres-1]->sicho-res[nres-2]->sicho;
  d/=d.norm();
  d.x()+=0.5;
  d.y()+=0.1;
  d.z()-=0.1;
  res[nres]->sicho=res[nres-1]->sicho+d*4.5;
  d.x()-=0.2;
  d.y()-=0.1;
  d.z()+=0.1;
  res[nres+1]->sicho=res[nres-1]->sicho+d*9.0;

  return nres+2;
}

void getFixedFromPDB(char *refpdb, AA **res, int* resinx, int* valid) {
  FILE *fptr=fopen(refpdb,"r");
  if (fptr==0) {
    fprintf(stderr,"cannot open reference PDB file %s\n",refpdb);
    return;
  }

  PDBEntry pdbEntry;

  while (!feof(fptr)) {
    if (pdbEntry.read(fptr)>0) {
      int num=pdbEntry.residueNumber();
      int rinx=resinx[num+100];
      if (rinx>=0 && !valid[rinx]) {
	if (strcmp(res[rinx]->name(),pdbEntry.residueName()) &&
	    (strcmp(res[rinx]->name(),"HIS") || 
	     (strcmp(pdbEntry.residueName(),"HSD") &&
	      strcmp(pdbEntry.residueName(),"HSE")))) {
	  fprintf(stderr,"residue %d: names >%s< and >%s< do not match\n",
		  resinx[num+100],res[rinx]->name(),pdbEntry.residueName());
	  exit(1);
	}
	res[rinx]->set(pdbEntry.atomName(),pdbEntry.coordinates());
      }
    }
  }

  fclose(fptr);
}

void readLibrary(Rotamer& rot, char *datafile) {
  ReducedLibRecord red;
  FILE *libptr=fopen(datafile,"rb");
  if (libptr==0) {
    fprintf(stderr,"cannot open library\n");
    exit(1);
  }
  int i,j,k;
  for (i=0; i<18; i++) {
    for (j=0; j<36*36; j++) {
      for (k=0; k<nrotamers[i]; k++) {
	if (!fread(&red,sizeof(ReducedLibRecord),1,libptr)) {
	  fprintf(stderr,"premature end of library file, aa:%d\n",red.aatype);
	  exit(1);
	}
	if (red.aatype!=i || red.phipsi!=j) {
	  fprintf(stderr,"invalid record in library: %d %d %lf %lf %lf %lf %lf\n",
		  red.aatype,red.phipsi,red.p,red.chi1,red.chi2,red.chi3,red.chi4);
	  exit(1);
	}
	rot.setRotamer(i,j,k,red);
      }
    }
  }
  fclose(libptr);
}

class EstCA {
public:
  Vector ca,mca,pca;
};

void getCAEstimates(char *cadname, AA **res, int nres, int* valid) {

  EstCA *est=new EstCA[20];

  double* ncnt=new double[nres];
  int i;
  for (i=0; i<nres; ncnt[i++]=0);

  FILE *fptr=fopen(cadname,"r");
  if (fptr==0) {
    fprintf(stderr,"cannot open CA data file >%s<\n",cadname);
    exit(1);
  }

  char line[1024];
  while(!feof(fptr)) {
    if (fgets(line,1023,fptr)!=0) {
      Field f(line," \n");
      int inx=atoi(f[0]);
      if (inx>=0 && inx<20) {
	est[inx].ca=Vector(atof(f[1]), atof(f[2]), atof(f[3]));
	est[inx].mca=Vector(atof(f[4]), atof(f[5]), atof(f[6]));
	est[inx].pca=Vector(atof(f[7]), atof(f[8]), atof(f[9]));
      } else {
	fprintf(stderr,"invalid data in CA data file: %s\n",line);
      }
    }
  }
  fclose(fptr);

  for (i=1; i<nres-1; i++) {
    Vector d12,d13,d23;
    
    d12=res[i]->sicho-res[i-1]->sicho;
    d12/=d12.norm();

    d23=res[i+1]->sicho-res[i]->sicho;
    d23/=d23.norm();

    d13=res[i+1]->sicho-res[i-1]->sicho;
    d13/=d13.norm();

    Vector ax,ay,az;

    ax=d13;

    az=d12.cross(d23);

    if (az.norm()<0.0001) {
      if (i+2<nres) {
	Vector d24=res[i+2]->sicho-res[i]->sicho;
	d24/=d24.norm();
	az=d12.cross(d24);
      } else {
	Vector d03=res[i+1]->sicho-res[i-2]->sicho;
	d03/=d03.norm();
	az=d12.cross(d03);
      }
    }
    az/=az.norm();

    ay=az.cross(ax);
    ay/=ay.norm();

    int inx=res[i]->numType;

    if (valid[i-1] && !res[i-1]->modelca) {
      res[i-1]->ca+=0.2*est[inx].mca.transformB(ax,ay,az);
      ncnt[i-1]+=0.2;
    }
    if (valid[i] && !res[i]->modelca) {
      res[i]->ca+=0.6*est[inx].ca.transformB(ax,ay,az);
      ncnt[i]+=0.6;
    }
    if (valid[i+1] && !res[i]->modelca) {
      res[i+1]->ca+=0.2*est[inx].pca.transformB(ax,ay,az);
      ncnt[i+1]+=0.2;
    }
  } 
  
  for (i=0; i<nres; i++) {
    if (valid[i] && !res[i]->modelca) {
      res[i]->ca/=ncnt[i];
      res[i]->ca+=res[i]->sicho;
    }
  }

  delete ncnt;
  delete[] est;
}

class EstCON {
public:
  int num;
  Vector c,n,o;
};

void getCONEstimates(char *pepname, AA **res, int nres, int* valid) {
  EstCON *est=new EstCON[20*16*16*32*2];

  FILE *fptr=fopen(pepname,"r");
  if (fptr==0) {
    fprintf(stderr,"cannot open CON data file\n");
    exit(1);
  }

  int i;
  for (i=0; i<20*16*16*32*2; i++) 
    est[i].num=0;
  
  int cnt=0;
  while(!feof(fptr)) {
    int inx,n;
    double xn,yn,zn,xc,yc,zc,xo,yo,zo;
    
    if (fscanf(fptr,"%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf",
	       &inx,&n,&xn,&yn,&zn,&xc,&yc,&zc,&xo,&yo,&zo)!=0) {
      if (inx>=0 && inx<20*16*16*32*2) {
	est[inx].num=n;
	est[inx].n=Vector(xn,yn,zn);
	est[inx].c=Vector(xc,yc,zc);
	est[inx].o=Vector(xo,yo,zo);
	cnt++;
      } 
    }
  }
  fclose(fptr);

  for (i=1; i<=nres-3; i++) {
    Vector d12,d13,d03,d23,d01,d02,t;

    d12=res[i+1]->ca-res[i]->ca;
    d13=res[i+2]->ca-res[i]->ca;
    d03=res[i+2]->ca-res[i-1]->ca;
    d02=res[i+1]->ca-res[i-1]->ca;
    d01=res[i]->ca-res[i-1]->ca;
    d23=res[i+2]->ca-res[i+1]->ca;
    
    Vector ax,ay,az;
    
    ax=d12.cross(d13);
    ax/=ax.norm();
    
    ay=d13.cross(ax);
    ay/=ay.norm();
    
    az=ax.cross(ay);
    ax/=ax.norm();

    int n1inx=vround((d02.norm()-4.0)/0.25);
    int n2inx=vround((d13.norm()-4.0)/0.25);
    int n3inx=vround((d03.norm()-4.0)/0.25);

    if (n1inx<0) n1inx=0;
    if (n2inx<0) n2inx=0;
    if (n3inx<0) n3inx=0;
    
    if (n1inx>=16) n1inx=15;
    if (n2inx>=16) n2inx=15;
    if (n3inx>=32) n3inx=31;

    t=d01.cross(d12);
    double v=t*d23;

    int cn3inx=(v<0)?n3inx+32:n3inx;

    int ninx=cn3inx*16*16+n2inx*16+n1inx;
    int inx=ninx*20+res[i]->numType;

    Vector ec,en,eo;
    if (est[inx].num<3) {
      int i1,i2,i3;

      int mini1,mini2,mini3;
      int maxi1,maxi2,maxi3;
      
      int add=1;
      int cnt=0;
      do {
	mini1=n1inx-add;
	mini2=n2inx-add;
	mini3=n3inx-add;
      
        if (mini1<0)  mini1=0;
	if (mini2<0)  mini2=0;
	if (mini3<0)  mini3=0;

	maxi1=n1inx+add;
	maxi2=n2inx+add;
	maxi3=n3inx+add;
	
	if (maxi1>=16) maxi1=15;
	if (maxi2>=16) maxi2=15;
	if (maxi3>=32) maxi3=31;

	for (i1=mini1; i1<=maxi1; i1++) {
	  for (i2=mini2; i2<=maxi2; i2++) {
	    for (i3=mini3; i3<=maxi3; i3++) {
	      int cn3=(v<0)?i3+32:i3;
	      ninx=cn3*16*16+i2*16+i1;
	      inx=ninx*20+res[i]->numType;
	      if (est[inx].num>0) {
		ec+=est[inx].c;
		eo+=est[inx].o;
		en+=est[inx].n;
		cnt++;
	      }
	    }
	  }
	}
	add++;
      } while (cnt==0);
      ec/=(double)cnt;
      en/=(double)cnt;
      eo/=(double)cnt;
    } else {
      ec=est[inx].c;
      eo=est[inx].o;
      en=est[inx].n;
    }

    if (valid[i]) {
      res[i]->c=res[i]->ca+ec.transformB(ax,ay,az);
      res[i]->n=res[i]->ca+en.transformB(ax,ay,az);
      res[i]->o=res[i]->ca+eo.transformB(ax,ay,az);
    }
  } 

  delete[] est;
}

double kinit=0.11;
//double kside=0.54;
double kside=1.0;
double kcaca=1.1;
double kcaxca=10.0;
double kpmf=0.25;

double kcond=5.0;
double kcond3=0.0;

double *grid02;
double *grid13;

void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	    double (*func)(double []), 
	    void (*dfunc)(double [], double []));

double casidemin[20] = { 3.250, 1.987, 1.984, 1.378, 2.210,
                         2.240, 2.699, 1.660, 1.930, 2.600,
                         2.120, 2.984, 1.405, 1.265, 1.455,
                         3.450, 3.350, 1.471, 0.763, 0.000 };

double casidemax[20] = { 4.100, 1.987, 1.984, 1.378, 2.785,
                         2.785, 2.699, 1.920, 2.097, 3.175,
                         2.640, 2.984, 1.405, 1.265, 1.455,
                         3.825, 3.350, 1.471, 0.763, 0.000 };

double potentialCA(double pca[]) {
  double e=0.0;
  double einit=0.0;
  double ecaca=0.0;
  double eside=0.0;
  double ecaxca=0.0;
  double epmf02=0.0;
  double epmf13=0.0;

  if (fixCA) return e;

  for (int i=1; i<=nres-3; i++) {
    Vector v,d,d01,d02;
    double dd,d0,te;

    v=Vector(pca[i*3+1],pca[i*3+2],pca[i*3+3]);

    // constraint around initial estimate

    d=v-res[i]->ca;
    dd=d*d;
    te=kinit*dd;

    e+=te;
    einit+=te;
      
    // constraint of CA - sidechain center of mass distance

    //    dd=v.norm();
    d=v-res[i]->sicho;
    dd=d.norm();
    
    d0=(dd<casidemin[res[i]->numType]) ? 
      dd-casidemin[res[i]->numType] : (dd>casidemax[res[i]->numType]) ? 
      dd-casidemax[res[i]->numType] : 0.0;

    te=kside*d0*d0;
    
    e+=te;    
    eside+=te;
      
    // constraint of CA-CA distance

    if (i>1) {
      Vector vm=Vector(pca[i*3-2],pca[i*3-1],pca[i*3]);
      
      //      v+=res[i]->sicho;
      //      vm+=res[i-1]->sicho;

      d01=v-vm;
      dd=d01.norm();

      if (res[i]->numType==12) {
	d0 = (dd<2.95)? dd-2.95 : ((dd>3.808)? dd-3.808 : 0.0);
      } else {
	d0=dd-3.808;
      }

      te=kcaca*d0*d0;
      
      e+=te;
      ecaca+=te;

      if (i<nres-3) {
	// constraint of CA - CA - CA angle through distance
	// between outer atoms
      
	Vector vp=Vector(pca[i*3+4],pca[i*3+5],pca[i*3+6]);
	
	d02=vp-vm;
	dd=d02.norm();

	double vmin=5.4;
        if (res[i]->numType==12) { //PRO
           vmin=4.5;
	}            
	
	d0=(dd<vmin) ? dd-vmin : (dd>7.3) ? dd-7.3 : 0.0;
	
	te=kcaxca*d0*d0;
	e+=te;
	ecaxca+=te;

	if (i<nres-4) {
	  Vector vpp=Vector(pca[i*3+7],pca[i*3+8],pca[i*3+9]);

	  Vector d12=vp-v;
	  Vector d13=vpp-v;
	  Vector d03=vpp-vm;
	  Vector d23=vpp-vp;

	  int i1=vround((d02.norm()-4.0)/0.2);
	  int i2=vround((d13.norm()-4.0)/0.2);

	  Vector t=d01.cross(d12);
	  double f=t*d23;

	  int i3=vround((d03.norm()-3.5)/0.2);
	  if (f>0) {
	    i3+=50;
	  }

	  if (i1>=0 && i1<20 && i3>=0 && i3<100) {
	    te=grid02[i1*100+i3];
	  } else {
	    te=6.0;
	  }
	  te*=kpmf;

	  e+=te;
	  epmf02+=te;

	  if (i2>=0 && i2<20 && i3>=0 && i3<100) {
	    te=grid13[i2*100+i3];
	  } else {
	    te=6.0;
	  }
	  te*=kpmf;

	  e+=te;
	  epmf13+=te;
	}
      }
    }
  }

  return e;
}

void gradientCA(double pca[], double xi[]) {
  int i;
  for (i=0; i<nres; i++) {
    xi[i*3+1]=0.0;
    xi[i*3+2]=0.0;
    xi[i*3+3]=0.0;
  }

  if (fixCA) return;
 
  Vector v,d,d01,d02;
  double dd,d0;

  for (i=1; i<=nres-3; i++) {
    v=Vector(pca[i*3+1],pca[i*3+2],pca[i*3+3]);
    
    // weak constraint around initial estimate

    d=v-res[i]->ca;
    dd=d*d;

    if (valid[i]) {
      xi[i*3+1]+=2.0*kinit*d.x();
      xi[i*3+2]+=2.0*kinit*d.y();
      xi[i*3+3]+=2.0*kinit*d.z();
    }
    
    // constraint of CA - sidechain center of mass distance
    
    d=v-res[i]->sicho;
    dd=d.norm();

    d0=(dd<casidemin[res[i]->numType]) ? 
      dd-casidemin[res[i]->numType] : (dd>casidemax[res[i]->numType]) ? 
      dd-casidemax[res[i]->numType] : 0.0;

    if (dd>0) {
      if (valid[i]) {
	xi[i*3+1]-=2.0*kside*d0*d.x()/dd;
	xi[i*3+2]-=2.0*kside*d0*d.y()/dd;
	xi[i*3+3]-=2.0*kside*d0*d.z()/dd;
      }
    }
  
    Vector vm;
    if (i>1) {
      // constraint of CA-CA distance
      vm=Vector(pca[i*3-2],pca[i*3-1],pca[i*3]);

      d01=v-vm;
      dd=d01.norm();

      if (res[i]->numType==12) {  // PRO
	d0 = (dd<2.95)? dd-2.95 : ((dd>3.808)? dd-3.808 : 0.0);
      } else {
	d0=dd-3.808;
      }
      
      if (dd>0) {
	if (valid[i]) {
	  xi[i*3+1]-=2.0*kcaca*d0*d01.x()/dd;
	  xi[i*3+2]-=2.0*kcaca*d0*d01.y()/dd;
	  xi[i*3+3]-=2.0*kcaca*d0*d01.z()/dd;
	}
	if (valid[i-1]) {
	  xi[(i-1)*3+1]+=2.0*kcaca*d0*d01.x()/dd;
	  xi[(i-1)*3+2]+=2.0*kcaca*d0*d01.y()/dd;
	  xi[(i-1)*3+3]+=2.0*kcaca*d0*d01.z()/dd;
	}
      }

      if (i<nres-3) {
	// constraint of CA - CA - CA angle through distance
	// between outer atoms

	Vector vp=Vector(pca[i*3+4],pca[i*3+5],pca[i*3+6]);
	
	d02=vp-vm;
	dd=d02.norm();

	double vmin=5.4;
        if (res[i]->numType==12) {  // PRO
           vmin=4.5;
	}            

	d0=(dd<vmin) ? dd-vmin : (dd>7.3) ? dd-7.3 : 0.0;
	
	if (dd>0.0) {
	  if (valid[i+1]) {
	    xi[(i+1)*3+1]-=2.0*kcaxca*d0*d02.x()/dd;
	    xi[(i+1)*3+2]-=2.0*kcaxca*d0*d02.y()/dd;
	    xi[(i+1)*3+3]-=2.0*kcaxca*d0*d02.z()/dd;
	  }
	  if (valid[i-1]) {
	    xi[(i-1)*3+1]+=2.0*kcaxca*d0*d02.x()/dd;
	    xi[(i-1)*3+2]+=2.0*kcaxca*d0*d02.y()/dd;
	    xi[(i-1)*3+3]+=2.0*kcaxca*d0*d02.z()/dd;
	  }
	}

	if (i<nres-4) {
	  Vector vpp=Vector(pca[i*3+7],pca[i*3+8],pca[i*3+9]);
	  
	  Vector d12=vp-v;
	  Vector d13=vpp-v;
	  Vector d03=vpp-vm;
	  Vector d23=vpp-vp;

	  int i1=vround((d02.norm()-4.0)/0.2);
	  int i2=vround((d13.norm()-4.0)/0.2);

	  Vector t=d01.cross(d12);
	  double f=t*d23;

	  int i3=vround((d03.norm()-3.5)/0.2);
	  if (f>0) 
	    i3+=50;

	  double te;
	  if (i1>0 && i1<19 && i3>0 && i3<99) {
	    te=grid02[i1*100+i3];
	    double t1m=grid02[(i1-1)*100+i3];
	    double t1p=grid02[(i1+1)*100+i3];
	    double t3m=grid02[i1*100+i3-1];
	    double t3p=grid02[i1*100+i3+1];

	    if (t1m<t1p && t1m<te) {
	      d0=te-t1m;
	    } else if (t1p<t1m && t1p<te) {
	      d0=t1p-te;
	    } else {
	      d0=0.0;
	    }

	    if (valid[i+1]) {	    
	      xi[(i+1)*3+1]-=kpmf*d0*d02.x();
	      xi[(i+1)*3+2]-=kpmf*d0*d02.y();
	      xi[(i+1)*3+3]-=kpmf*d0*d02.z();
	    }
	  
	    if (valid[i-1]) {	    
	      xi[(i-1)*3+1]+=kpmf*d0*d02.x();
	      xi[(i-1)*3+2]+=kpmf*d0*d02.y();
	      xi[(i-1)*3+3]+=kpmf*d0*d02.z();
	    }

	    if (t3m<t3p && t3m<te) {
	      d0=te-t3m;
	    } else if (t3p<t3m && t3p<te) {
	      d0=t3p-te;
	    } else {
	      d0=0.0;
	    }
	    
	    if (valid[i+2]) {	    
	      xi[(i+2)*3+1]-=kpmf*d0*d03.x();
	      xi[(i+2)*3+2]-=kpmf*d0*d03.y();
	      xi[(i+2)*3+3]-=kpmf*d0*d03.z();
	    }
	  
	    if (valid[i-1]) {	    
	      xi[(i-1)*3+1]+=kpmf*d0*d03.x();
	      xi[(i-1)*3+2]+=kpmf*d0*d03.y();
	      xi[(i-1)*3+3]+=kpmf*d0*d03.z();
	    }
	  }

	  if (i2>0 && i2<19 && i3>0 && i3<99) {
	    te=grid13[i2*100+i3];
	    double t2m=grid13[(i2-1)*100+i3];
	    double t2p=grid13[(i2+1)*100+i3];
	    double t3m=grid13[i2*100+i3-1];
	    double t3p=grid13[i2*100+i3+1];

	    if (t2m<t2p && t2m<te) {
	      d0=te-t2m;
	    } else if (t2p<t2m && t2p<te) {
	      d0=t2p-te;
	    } else {
	      d0=0.0;
	    }
	    
	    if (valid[i+2]) {	    
	      xi[(i+2)*3+1]-=kpmf*d0*d13.x();
	      xi[(i+2)*3+2]-=kpmf*d0*d13.y();
	      xi[(i+2)*3+3]-=kpmf*d0*d13.z();
	    }
	  
	    if (valid[i]) {	    
	      xi[i*3+1]+=kpmf*d0*d13.x();
	      xi[i*3+2]+=kpmf*d0*d13.y();
	      xi[i*3+3]+=kpmf*d0*d13.z();
	    }

	    if (t3m<t3p && t3m<te) {
	      d0=te-t3m;
	    } else if (t3p<t3m && t3p<te) {
	      d0=t3p-te;
	    } else {
	      d0=0.0;
	    }
	    
	    if (valid[i+2]) {	    
	      xi[(i+2)*3+1]-=kpmf*d0*d03.x();
	      xi[(i+2)*3+2]-=kpmf*d0*d03.y();
	      xi[(i+2)*3+3]-=kpmf*d0*d03.z();
	    }
	  
	    if (valid[i-1]) {	    
	      xi[(i-1)*3+1]+=kpmf*d0*d03.x();
	      xi[(i-1)*3+2]+=kpmf*d0*d03.y();
	      xi[(i-1)*3+3]+=kpmf*d0*d03.z();
	    }
	  }
	}
      }
    }
  }
}

double potentialCON(double pb[]) {
  double e=potentialCA(pb);
  double econd=0.0;
  double econd3=0.0;

  Vector d;
  double d0,te;
  
  for (int i=2; i<=nres-3; i++) {
    Vector vca=Vector(pb[i*3+1],pb[i*3+2],pb[i*3+3]);
    Vector vn=Vector(pb[nres*3+i*3+1],pb[nres*3+i*3+2],pb[nres*3+i*3+3]);
    Vector vc=Vector(pb[nres*3*2+i*3+1],pb[nres*3*2+i*3+2],pb[nres*3*2+i*3+3]);
    Vector vo=Vector(pb[nres*3*3+i*3+1],pb[nres*3*3+i*3+2],pb[nres*3*3+i*3+3]);

    Vector vcm=Vector(pb[nres*3*2+(i-1)*3+1],pb[nres*3*2+(i-1)*3+2],
		      pb[nres*3*2+(i-1)*3+3]);
    Vector vom=Vector(pb[nres*3*3+(i-1)*3+1],pb[nres*3*3+(i-1)*3+2],
		      pb[nres*3*3+(i-1)*3+3]);
    Vector vcam=Vector(pb[(i-1)*3+1],pb[(i-1)*3+2],pb[(i-1)*3+3]);

    d=vc-vo;
    d0=d.norm()-1.2325;

    te=kcond*d0*d0;
    e+=te;
    econd+=te;

    d=vca-vn;
    d0=d.norm()-1.460;
    te=kcond*d0*d0;
    e+=te;
    econd+=te;

    d=vca-vc;
    d0=d.norm()-1.525;
    te=kcond*d0*d0;
    e+=te;
    econd+=te;
    
    d=vn-vcm;
    d0=d.norm()-1.33;
    te=kcond*d0*d0;
    e+=te;
    econd+=te;
    
    d=vc-vn;
    d0=d.norm()-2.46;
    te=kcond3*d0*d0;
    e+=te;
    econd3+=te;
    
    d=vca-vcm;
    d0=d.norm()-2.435;
    te=kcond3*d0*d0;
    e+=te;
    econd3+=te;

    d=vn-vcam;
    d0=d.norm()-2.4325;
    te=kcond3*d0*d0;
    e+=te;
    econd3+=te;

    d=vo-vca;
    d0=d.norm()-2.4;
    te=kcond3*d0*d0;
    e+=te;
    econd3+=te;

    d=vn-vom;
    d0=d.norm()-2.25;
    te=kcond3*d0*d0;
    e+=te;
    econd3+=te;
  }

  return e;
}

void gradientCON(double pb[], double xi[]) {
  gradientCA(pb,xi);

  int i,k;
  for (k=1; k<=3; k++) {
    for (i=0; i<nres; i++) {
      xi[nres*3*k+i*3+1]=0.0;
      xi[nres*3*k+i*3+2]=0.0;
      xi[nres*3*k+i*3+3]=0.0;
    }
  }

  Vector d;
  double dd,d0;

  for (i=2; i<=nres-3; i++) {
    Vector vca=Vector(pb[i*3+1],pb[i*3+2],pb[i*3+3]);
    Vector vn=Vector(pb[nres*3+i*3+1],pb[nres*3+i*3+2],pb[nres*3+i*3+3]);
    Vector vc=Vector(pb[nres*3*2+i*3+1],pb[nres*3*2+i*3+2],pb[nres*3*2+i*3+3]);
    Vector vo=Vector(pb[nres*3*3+i*3+1],pb[nres*3*3+i*3+2],pb[nres*3*3+i*3+3]);

    Vector vcm=Vector(pb[nres*3*2+(i-1)*3+1],pb[nres*3*2+(i-1)*3+2],
		      pb[nres*3*2+(i-1)*3+3]);
    Vector vom=Vector(pb[nres*3*3+(i-1)*3+1],pb[nres*3*3+(i-1)*3+2],
		      pb[nres*3*3+(i-1)*3+3]);
    Vector vcam=Vector(pb[(i-1)*3+1],pb[(i-1)*3+2],pb[(i-1)*3+3]);

    d=vc-vo;
    dd=d.norm();
    d0=dd-1.2325;

    if (dd>0) {
      if (valid[i]) {	    
	xi[nres*3*2+i*3+1]-=2.0*kcond*d0*d.x()/dd;
	xi[nres*3*2+i*3+2]-=2.0*kcond*d0*d.y()/dd;
	xi[nres*3*2+i*3+3]-=2.0*kcond*d0*d.z()/dd;

	xi[nres*3*3+i*3+1]+=2.0*kcond*d0*d.x()/dd;
	xi[nres*3*3+i*3+2]+=2.0*kcond*d0*d.y()/dd;
	xi[nres*3*3+i*3+3]+=2.0*kcond*d0*d.z()/dd;
      }
    }

    d=vca-vn;
    dd=d.norm();
    d0=dd-1.460;


    if (dd>0) {
      if (valid[i]) {	    
	if (!fixCA) {
	  xi[i*3+1]-=2.0*kcond*d0*d.x()/dd;
	  xi[i*3+2]-=2.0*kcond*d0*d.y()/dd;
	  xi[i*3+3]-=2.0*kcond*d0*d.z()/dd;
	}

	xi[nres*3*1+i*3+1]+=2.0*kcond*d0*d.x()/dd;
	xi[nres*3*1+i*3+2]+=2.0*kcond*d0*d.y()/dd;
	xi[nres*3*1+i*3+3]+=2.0*kcond*d0*d.z()/dd;
      }
    }

    d=vca-vc;
    dd=d.norm();
    d0=dd-1.525;


    if (dd>0) {
      if (valid[i]) {
	if (!fixCA) {
	  xi[i*3+1]-=2.0*kcond*d0*d.x()/dd;
	  xi[i*3+2]-=2.0*kcond*d0*d.y()/dd;
	  xi[i*3+3]-=2.0*kcond*d0*d.z()/dd;
	}

	xi[nres*3*2+i*3+1]+=2.0*kcond*d0*d.x()/dd;
	xi[nres*3*2+i*3+2]+=2.0*kcond*d0*d.y()/dd;
	xi[nres*3*2+i*3+3]+=2.0*kcond*d0*d.z()/dd;
      }
    }

    d=vn-vcm;
    dd=d.norm();
    d0=dd-1.33;

    if (dd>0) {
      if (valid[i]) {	    
	xi[nres*3*1+i*3+1]-=2.0*kcond*d0*d.x()/dd;
	xi[nres*3*1+i*3+2]-=2.0*kcond*d0*d.y()/dd;
	xi[nres*3*1+i*3+3]-=2.0*kcond*d0*d.z()/dd;
      }

      if (valid[i-1]) {	    
	xi[nres*3*2+(i-1)*3+1]+=2.0*kcond*d0*d.x()/dd;
	xi[nres*3*2+(i-1)*3+2]+=2.0*kcond*d0*d.y()/dd;
	xi[nres*3*2+(i-1)*3+3]+=2.0*kcond*d0*d.z()/dd;
      }
    }
    
    d=vc-vn;
    dd=d.norm();
    d0=dd-2.46;

    if (dd>0) {
      if (valid[i]) {	    
	xi[nres*3*2+i*3+1]-=2.0*kcond3*d0*d.x()/dd;
	xi[nres*3*2+i*3+2]-=2.0*kcond3*d0*d.y()/dd;
	xi[nres*3*2+i*3+3]-=2.0*kcond3*d0*d.z()/dd;

	xi[nres*3*1+i*3+1]+=2.0*kcond3*d0*d.x()/dd;
	xi[nres*3*1+i*3+2]+=2.0*kcond3*d0*d.y()/dd;
	xi[nres*3*1+i*3+3]+=2.0*kcond3*d0*d.z()/dd;
      }
    }

    d=vca-vcm;
    dd=d.norm();
    d0=dd-2.435;

    if (dd>0) {
      if (valid[i]) {	    
	if (!fixCA) {
	  xi[i*3+1]-=2.0*kcond3*d0*d.x()/dd;
	  xi[i*3+2]-=2.0*kcond3*d0*d.y()/dd;
	  xi[i*3+3]-=2.0*kcond3*d0*d.z()/dd;
	}
      }

      if (valid[i-1]) {	    
	xi[nres*3*2+(i-1)*3+1]+=2.0*kcond3*d0*d.x()/dd;
	xi[nres*3*2+(i-1)*3+2]+=2.0*kcond3*d0*d.y()/dd;
	xi[nres*3*2+(i-1)*3+3]+=2.0*kcond3*d0*d.z()/dd;
      }
    }

    d=vn-vcam;
    dd=d.norm();
    d0=dd-2.4325;

    if (dd>0) {
      if (valid[i]) {	    
	xi[nres*3*1+i*3+1]-=2.0*kcond3*d0*d.x()/dd;
	xi[nres*3*1+i*3+2]-=2.0*kcond3*d0*d.y()/dd;
	xi[nres*3*1+i*3+3]-=2.0*kcond3*d0*d.z()/dd;
      }

      if (valid[i-1]) {	    
	if (!fixCA) {
	  xi[(i-1)*3+1]+=2.0*kcond3*d0*d.x()/dd;
	  xi[(i-1)*3+2]+=2.0*kcond3*d0*d.y()/dd;
	  xi[(i-1)*3+3]+=2.0*kcond3*d0*d.z()/dd;
	}
      }
    }

    d=vo-vca;
    dd=d.norm();
    d0=dd-2.4;
    
    if (dd>0) {
      if (valid[i]) {	    
	xi[nres*3*3+i*3+1]-=2.0*kcond3*d0*d.x()/dd;
	xi[nres*3*3+i*3+2]-=2.0*kcond3*d0*d.y()/dd;
	xi[nres*3*3+i*3+3]-=2.0*kcond3*d0*d.z()/dd;

	if (!fixCA) {
	  xi[i*3+1]+=2.0*kcond3*d0*d.x()/dd;
	  xi[i*3+2]+=2.0*kcond3*d0*d.y()/dd;
	  xi[i*3+3]+=2.0*kcond3*d0*d.z()/dd;
	}
      }
    }

    d=vn-vom;
    dd=d.norm();
    d0=dd-2.25;

    if (dd>0) {
      if (valid[i]) {	    
	xi[nres*3*1+i*3+1]-=2.0*kcond3*d0*d.x()/dd;
	xi[nres*3*1+i*3+2]-=2.0*kcond3*d0*d.y()/dd;
	xi[nres*3*1+i*3+3]-=2.0*kcond3*d0*d.z()/dd;
      }

      if (valid[i-1]) {	    
	xi[nres*3*3+(i-1)*3+1]+=2.0*kcond3*d0*d.x()/dd;
	xi[nres*3*3+(i-1)*3+2]+=2.0*kcond3*d0*d.y()/dd;
	xi[nres*3*3+(i-1)*3+3]+=2.0*kcond3*d0*d.z()/dd;
      }
    }
  }
}


void minimizeCA(AA **res, int nres, int limited) {
  int i;
  double *pca=new double[nres*3];
  for (i=0; i<nres; i++) {
    pca[i*3]=res[i]->ca.x();
    pca[i*3+1]=res[i]->ca.y();
    pca[i*3+2]=res[i]->ca.z();
  }

  int iter; 
  double fret;

  frprmn(pca-1,nres*3,0.1,&iter,&fret,
  	 potentialCA,gradientCA);

  kcaca=200.0;
  kcaxca=200.0;

  frprmn(pca-1,nres*3,0.0000001,&iter,&fret,
  	 potentialCA,gradientCA);
  
  Vector ca(0.0,0.0,0.0);
  
  for (i=0; i<nres; i++) {
    res[i]->ca=Vector(pca[i*3],pca[i*3+1],pca[i*3+2]);
    ca+=res[i]->ca-res[i]->sicho;
  }

  if (!limited) {
    ca/=nres;
  
    for (i=0; i<nres; i++) {
      res[i]->ca+=ca;
    }
  }

  delete pca;
}

void minimizeCON(AA **res, int nres) {
  int i;
  double *pb=new double[nres*3*4];

  for (i=0; i<nres; i++) {
    pb[i*3]=res[i]->ca.x();
    pb[i*3+1]=res[i]->ca.y();
    pb[i*3+2]=res[i]->ca.z();

    pb[nres*3+i*3]=res[i]->n.x();
    pb[nres*3+i*3+1]=res[i]->n.y();
    pb[nres*3+i*3+2]=res[i]->n.z();

    pb[nres*3*2+i*3]=res[i]->c.x();
    pb[nres*3*2+i*3+1]=res[i]->c.y();
    pb[nres*3*2+i*3+2]=res[i]->c.z();

    pb[nres*3*3+i*3]=res[i]->o.x();
    pb[nres*3*3+i*3+1]=res[i]->o.y();
    pb[nres*3*3+i*3+2]=res[i]->o.z();
  }

  int iter; 
  double fret;

  frprmn(pb-1,nres*3*4,0.0000001,&iter,&fret,
    	 potentialCON,gradientCON);

  for (i=0; i<nres; i++) {
    res[i]->ca=Vector(pb[i*3],pb[i*3+1],pb[i*3+2]);
    res[i]->n=Vector(pb[nres*3+i*3],pb[nres*3+i*3+1],pb[nres*3+i*3+2]);
    res[i]->c=Vector(pb[nres*3*2+i*3],pb[nres*3*2+i*3+1],pb[nres*3*2+i*3+2]);    
    res[i]->o=Vector(pb[nres*3*3+i*3],pb[nres*3*3+i*3+1],pb[nres*3*3+i*3+2]);
  }

  delete pb;
}

void printBackbone(AA **res, int nres) {  
  int iatom=1;
  for (int i=1; i<nres-2; i++) {
    res[i]->printBack(iatom);
  }
}

void printComplete(AA **res, int nres) {
  int natom=1;
  for (int i=1; i<nres-2; i++) {
    res[i]->print(natom);
  }
}

void readPMFGrid(char *fname, double*& g) {
  g=new double[20*100];

  FILE *fptr=fopen(fname,"r");
  if (fptr==0) {
    fprintf(stderr,"cannot open grid file %s\n",fname);
    exit(1);
  }
  
  int i,j;
  double v;

  for (i=0; i<20*100; i++) {
    g[i]=6.0;
  }

  while(!feof(fptr)) {
    if (fscanf(fptr,"%d%d%lf",&i,&j,&v)) {
      if (i>=0 && i<20 && j>=0 && j<100) 
	g[i*100+j]=v;
    }
  }
  fclose(fptr);
}


void placeSideChains(AA **res, int nres, Rotamer& rot, int *valid) {
  int i,j,k;
  for (i=1; i<nres-2; i++) {
    if (valid[i]) {
      double minx[4];

      int aainx=res[i]->numType;
      if (aainx>=18) {
	res[i]->setChi(0,0,0,0);
      } else {
	double phi,psi;
    
	if (i>1) 
	  phi=dihedral(res[i-1]->c,res[i]->n,res[i]->ca,res[i]->c);
	else
	  phi=-90;
      
	if (i<nres) 
	  psi=dihedral(res[i]->n,res[i]->ca,res[i]->c,res[i+1]->n);
	else
	  psi=0.0;
      
	int phipsi=getPhiPsiIndex(phi,psi);
	AA *a=getAA(res[i]->name());
	*a=*(res[i]);

	double x[4];

	double mind=1.0E99;
	for (j=0; j<nrotamers[aainx]; j++) {
	  x[0]=rot.chi1(aainx,phipsi,j);
	  x[1]=rot.chi2(aainx,phipsi,j);
	  x[2]=rot.chi3(aainx,phipsi,j);
	  x[3]=rot.chi4(aainx,phipsi,j);

	  a->setChi(x[0],x[1],x[2],x[3]);
	  a->makeSide();
	  Vector tcofm=a->cofm();
	  Vector delta=tcofm-res[i]->sicho;
	  double d=delta.norm();
   
	  if (d<mind) {
	    mind=d;
	    for (k=0; k<4; k++) 
	      minx[k]=x[k];
	  } 
	}
	
	double add=2.0;
	
	for (k=0; k<4; k++) 
	  x[k]=minx[k];

	for (j=0; j<a->nx; j++) {
	  double mindp,mindn;
	  double minxp,minxn;
	  double d;

	  d=mind;
	  int cnt=0;
	  do {
	    mindp=d;
	    minxp=x[j];
	    x[j]+=add;
	    a->setChi(x[0],x[1],x[2],x[3]);
	    a->makeSide();
	    Vector tcofm=a->cofm();
	    Vector delta=tcofm-res[i]->sicho;
	    d=delta.norm();
	  } while(d<mindp && ++cnt<10);
	  
	  x[j]=minx[j];
	  d=mind;
	  cnt=0;
	  do {
	    mindn=d;
	    minxn=x[j];
	    x[j]-=add;
	    a->setChi(x[0],x[1],x[2],x[3]);
	    a->makeSide();
	    Vector tcofm=a->cofm();
	    Vector delta=tcofm-res[i]->sicho;
	    d=delta.norm();
	  } while(d<mindn && ++cnt<10);

	  x[j] = minx[j] = (mindn<mindp) ? minxn : minxp;
	  mind = (mindn<mindp) ? mindn : mindp;
	}
        delete a;
      }
      res[i]->setChi(minx[0],minx[1],minx[2],minx[3]);
      res[i]->makeSide();
    }
  }
}

int main(int argc, char **argv) {
  char refpdb[256];
  *refpdb=0;

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

  int *resinx=new int[MAXRES];
  for (i=0; i<MAXRES; resinx[i++]=-1);
  for (i=0; i<MAXRES; valid[i++]=1);

  int limmin[MAXRES];
  int limmax[MAXRES];
  int nlim=0;

  int backboneOnly=0;
  fixCA=0;

  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i],"-l")) {

      Field f(argv[++i],"=");
      for (j=0; j<f.number(); j++) {
	Field g(f[j],":");
	limmin[nlim]=atoi(g[0]);
	if (g.number()>1) {
	  limmax[nlim++]=atoi(g[1]);
	} else {
	  limmax[nlim]=limmin[nlim];
	  nlim++;
	}
      }

      strcpy(refpdb,argv[++i]);
    } else if (!strcmp(argv[i],"-d")) {
      strcpy(datadir,argv[++i]);
    } else if (!strcmp(argv[i],"-backonly")) {
      backboneOnly=1;
    } else if (!strcmp(argv[i],"-fixca")) {
      fixCA=1;
    } else if (!strcmp(argv[i],"-help")) {
      usage();
    } else {
      strcpy(inputfile,argv[i]);
      break;
    }
  }

  nres=readInputFile(inputfile,res,resinx);

  if (nlim>0) {
    for (i=0; i<MAXRES; valid[i++]=0);
    for (i=0; i<nlim; i++) 
      for (j=limmin[i]; j<=limmax[i]; j++) {
	if (resinx[j+100]>=0) {
	  valid[resinx[j+100]]=1;
	}
      }
    if (valid[1])
      valid[0]=1;
    if (nres>=3 && valid[nres-3])
      valid[nres-2]=valid[nres-1]=1;

    getFixedFromPDB(refpdb,res,resinx,valid);    
  }


  char fname[512];

  strcpy(fname,datadir);strcat(fname,"/");
  strcat(fname,"ca.lib");
  getCAEstimates(fname,res,nres,valid);

  strcpy(fname,datadir);strcat(fname,"/");
  strcat(fname,"grid02");
  readPMFGrid(fname,grid02);
  strcpy(fname,datadir);strcat(fname,"/");
  strcat(fname,"grid13");
  readPMFGrid(fname,grid13);

  if (!fixCA) 
    minimizeCA(res,nres,(nlim>0));

  strcpy(fname,datadir);strcat(fname,"/");
  strcat(fname,"con.lib");
  getCONEstimates(fname,res,nres,valid);

  minimizeCON(res,nres);

  if (backboneOnly)  
    printBackbone(res,nres);
  else {
    char scdatafile[512];
    strcpy(scdatafile,datadir);
    strcat(scdatafile,"/");
    strcat(scdatafile,SCWRLLIB);
    
    Rotamer rot(20,36*36,nrotamers);
    readLibrary(rot,scdatafile);

    placeSideChains(res,nres,rot,valid);

    printComplete(res,nres);
  }

  delete[] resinx;

  return 0;
}

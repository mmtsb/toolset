#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "vector.h"
#include "aa.h"

// *** AA ******

void AA::setChi(double chi1, double chi2, double chi3, double chi4) {
  x[0]=chi1;
  x[1]=chi2;
  x[2]=chi3;
  x[3]=chi4;
}

int AA::setBack(char *name, Vector v) {
  if (!strcmp(name,"CA")) 
    ca=v;
  else if (!strcmp(name,"C"))
    c=v;
  else if (!strcmp(name,"N"))
    n=v;
  else if (!strcmp(name,"O") || !strcmp(name,"OT1"))
    o=v;
  else 
    return -1;
  return 0;
}

/*
 * The following routines nmcrosprod, viewat are from the SCWRL code
 * by R. L. Dunbrack, RL_Dunbrack@fccc.edu
 *
 * They are originally from Paul Bash.  The math involved in
 * putting an arbitrary vector onto the z axis may be found in
 * Newman & Sproul, "Principles of Interactive Computer Graphics".
 *
 * Further modifications by Don Kneller 12/7/88) */

void AA::nmcrosprod(double  x1, double y1, double z1, 
		    double  x2, double y2, double z2, 
		    double& x3, double& y3, double& z3) {
  double        dis;		/* length of crossproduct vector r1 x r2 */
  double	x, y, z;
  
  x = y1*z2 - y2*z1;
  y = z1*x2 - z2*x1;
  z = x1*y2 - x2*y1;
  
  dis = sqrt((x*x + y*y + z*z));
  
  x3 = x / dis;
  y3 = y / dis;
  z3 = z / dis;
}


void AA::viewat(double **M, double **invM, 
		double *P1, double *P2, double *P3) {
  double	d12;
  double	P120, P121, P122, P130, P131, P132;
  
  P120 = P2[0] - P1[0];
  P121 = P2[1] - P1[1];
  P122 = P2[2] - P1[2];
  P130 = P3[0] - P1[0];
  P131 = P3[1] - P1[1];
  P132 = P3[2] - P1[2];
  d12 =sqrt(P120*P120 + P121*P121 + P122*P122);
  invM[2][0] = M[0][2] = P120/d12;
  invM[2][1] = M[1][2] = P121/d12;
  invM[2][2] = M[2][2] = P122/d12;
  
  nmcrosprod(P130, P131, P132, P120, P121, P122,
	     M[0][0], M[1][0], M[2][0]);
  invM[0][0] = M[0][0];
  invM[0][1] = M[1][0];
  invM[0][2] = M[2][0];
  
  M[0][3] = invM[0][3] = 0.;
  M[1][3] = invM[1][3] = 0.;
  M[2][3] = invM[2][3] = 0.;
  
  nmcrosprod(M[0][2], M[1][2], M[2][2], M[0][0], M[1][0], M[2][0],
	     M[0][1], M[1][1], M[2][1]);
  invM[1][0] = M[0][1];
  invM[1][1] = M[1][1];
  invM[1][2] = M[2][1];
  
  invM[3][0] = P1[0];
  invM[3][1] = P1[1];
  invM[3][2] = P1[2];
  invM[3][3] = 1.0;
  
  M[3][0] = -P1[0]*M[0][0] - P1[1]*M[1][0] - P1[2]*M[2][0];
  M[3][1] = -P1[0]*M[0][1] - P1[1]*M[1][1] - P1[2]*M[2][1];
  M[3][2] = -P1[0]*M[0][2] - P1[1]*M[1][2] - P1[2]*M[2][2];
  M[3][3] = 1.0;
}

void AA::mtxmult(double **mat, Vector& w) {
  int 		i, j;
  double 	oldvect[4];
  double	newvect[4];
  
  oldvect[0] = w.x();
  oldvect[1] = w.y();
  oldvect[2] = w.z();
  oldvect[3] = 1;
  
  for (i = 0; i < 4; i++) {
    newvect[i] = 0.0;
  }
  
  for (i = 0; i < 4; i++) {
    for(j = 0; j < 4; j++) {
      newvect[i] += oldvect[j] * mat[j][i];
    }
  }
  
  w=Vector(newvect[0],newvect[1],newvect[2]);
}  


void AA::build(Vector& w, double bond, double angle, double dih, 
	       Vector& v1, Vector &v2, Vector &v3) {
  angle*=M_PI/180.0;
  dih*=M_PI/180.0;

  w=Vector(bond*sin(angle)*sin(dih),bond*sin(angle)*cos(dih),bond*cos(angle));
 
  double *mat[4];
  double *invmat[4];
  
  int i;
  for (i=0; i<4; i++) {
    mat[i]=new double[4];
    invmat[i]=new double[4];
  }

  double p1[4],p2[4],p3[4];

  p1[0]=v1.x(); 
  p1[1]=v1.y();
  p1[2]=v1.z();

  p2[0]=v2.x(); 
  p2[1]=v2.y();
  p2[2]=v2.z();

  p3[0]=v3.x(); 
  p3[1]=v3.y();
  p3[2]=v3.z();
  
  p1[3]=p2[3]=p3[3]=1.0;

  viewat(mat,invmat,p1,p2,p3);

  mtxmult(invmat,w);

  for (i=0; i<4; i++) {
    delete mat[i];
    delete invmat[i];
  }
}

void AA::printPDBLine(Vector& v, char *aname, int natom, char *resname) {
  fprintf(stdout,
	  "ATOM %6d  %-3s %-4s%c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",
	  natom,aname,resname,pdbchain,pdbnum,v.x(),v.y(),v.z(),0.0,0.0,pdbsegid);
}

void AA::printBack(int& natom) {
  printPDBLine(n,"N",natom++,name());
  printPDBLine(ca,"CA",natom++,name());
  printPDBLine(c,"C",natom++,name());
  printPDBLine(o,"O",natom++,name());
}

void AA::addVectorBack(Vector *v, int& natom) {
  v[natom++]=n;
  v[natom++]=ca;
  v[natom++]=c;
  v[natom++]=o;
}

// *** Arg ******

int Arg::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG"))
      cg=v;
    else if (!strcmp(name,"CD"))
      cd=v;
    else if (!strcmp(name,"NE"))
      ne=v;
    else if (!strcmp(name,"CZ"))
      cz=v;
    else if (!strcmp(name,"NH1"))
      nh1=v;
    else if (!strcmp(name,"NH2"))
      nh2=v;
    else 
      return -1;
  }
  return 0;
}

Vector Arg::cofm() {
  Vector v;
  v=ca+cb+cg+cd+ne+cz+nh1+nh2;
  v/=8.0;
  return v;
}

void Arg::makeSide() {
  build(cb,  1.530, 110.5, -120.0, ca,  n,  c);
  build(cg,  1.520, 114.1,     x[0], cb, ca,  n);
  build(cd,  1.520, 111.3,     x[1], cg, cb, ca);
  build(ne,  1.461, 112.0,     x[2], cd, cg, cb);
  build(cz,  1.329, 124.2,     x[3], ne, cd, cg);
  build(nh1, 1.326, 120.0,  180.0, cz, ne, cd);
  build(nh2, 1.326, 120.0,    0.0, cz, ne, cd);
}

void Arg::print(int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg,"CG",natom++,name());
  printPDBLine(cd,"CD",natom++,name());
  printPDBLine(ne,"NE",natom++,name());
  printPDBLine(cz,"CZ",natom++,name());
  printPDBLine(nh1,"NH1",natom++,name());
  printPDBLine(nh2,"NH2",natom++,name());
}

void Arg::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg;
  v[natom++]=cd;
  v[natom++]=ne;
  v[natom++]=cz;
  v[natom++]=nh1;
  v[natom++]=nh2;
}


// *** Asn ******

int Asn::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG"))
      cg=v;
    else if (!strcmp(name,"OD1"))
      od1=v;
    else if (!strcmp(name,"ND2"))
      nd2=v;
    else 
      return -1;
  }
  return 0;
}

Vector Asn::cofm() {
  Vector v;
  v=ca+cb+cg+od1+nd2;
  v/=5.0;
  return v;
}

void Asn::makeSide() {
  build(cb,  1.530, 110.5, -120.0, ca,  n,  c);
  build(cg,  1.516, 112.6,       x[0], cb, ca,  n);
  build(od1, 1.231, 120.8,       x[1], cg, cb, ca);
  build(nd2, 1.328, 116.4, x[1]-180.0, cg, cb, ca);
}

void Asn::print(int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg,"CG",natom++,name());
  printPDBLine(od1,"OD1",natom++,name());
  printPDBLine(nd2,"ND2",natom++,name());
}

void Asn::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg;
  v[natom++]=od1;
  v[natom++]=nd2;
}


// *** Asp ******

int Asp::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG"))
      cg=v;
    else if (!strcmp(name,"OD1"))
      od1=v;
    else if (!strcmp(name,"OD2"))
      od2=v;
    else 
      return -1;
  }
  return 0;
}

Vector Asp::cofm() {
  Vector v;
  v=ca+cb+cg+od1+od2;
  v/=5.0;
  return v;
}

void Asp::makeSide() {
  build(cb,  1.530, 110.5,   -120.0, ca,  n,  c);
  build(cg,  1.516, 112.6,       x[0], cb, ca,  n);
  build(od1, 1.249, 118.4,       x[1], cg, cb, ca);
  build(od2, 1.249, 118.4, x[1]-180.0, cg, cb, ca);
}

void Asp::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg,"CG",natom++,name());
  printPDBLine(od1,"OD1",natom++,name());
  printPDBLine(od2,"OD2",natom++,name());
}

void Asp::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg;
  v[natom++]=od1;
  v[natom++]=od2;
}

// *** Cys ******

int Cys::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"SG"))
      sg=v;
    else 
      return -1;
  }
  return 0;
}

Vector Cys::cofm() {
  Vector v;
  v=ca+cb+sg;
  v/=3.0;
  return v;
}

void Cys::makeSide() {
  build(cb,  1.530, 110.5,   -120.0, ca,  n,  c);
  build(sg,  1.808, 114.4,       x[0], cb, ca,  n);
}

void Cys::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(sg,"SG",natom++,name());
}

void Cys::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=sg;
}


// *** Gln ******

int Gln::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG"))
      cg=v;
    else if (!strcmp(name,"CD"))
      cd=v;
    else if (!strcmp(name,"OE1"))
      oe1=v;
    else if (!strcmp(name,"NE2"))
      ne2=v;
    else 
      return -1;
  }
  return 0;
}

Vector Gln::cofm() {
  Vector v;
  v=ca+cb+cg+cd+oe1+ne2;
  v/=6.0;
  return v;
}

void Gln::makeSide() {
  build(cb,  1.530, 110.5,   -120.0, ca,  n,  c);
  build(cg,  1.520, 114.1,       x[0], cb, ca,  n);
  build(cd,  1.516, 112.6,       x[1], cg, cb, ca);
  build(oe1, 1.231, 120.8,       x[2], cd, cg, cb);
  build(ne2, 1.328, 116.4, x[2]+180.0, cd, cg, cb);
}

void Gln::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg,"CG",natom++,name());
  printPDBLine(cd,"CD",natom++,name());
  printPDBLine(oe1,"OE1",natom++,name());
  printPDBLine(ne2,"NE2",natom++,name());
}

void Gln::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg;
  v[natom++]=cd;
  v[natom++]=oe1;
  v[natom++]=ne2;
}

// *** Glu ******

int Glu::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG"))
      cg=v;
    else if (!strcmp(name,"CD"))
      cd=v;
    else if (!strcmp(name,"OE1"))
      oe1=v;
    else if (!strcmp(name,"OE2"))
      oe2=v;
    else 
      return -1;
  }
  return 0;
}

Vector Glu::cofm() {
  Vector v;
  v=ca+cb+cg+cd+oe1+oe2;
  v/=6.0;
  return v;
}

void Glu::makeSide() {
  build(cb,  1.530, 110.5,   -120.0, ca,  n,  c);
  build(cg,  1.520, 114.1,       x[0], cb, ca,  n);
  build(cd,  1.516, 112.6,       x[1], cg, cb, ca);
  build(oe1, 1.249, 118.4,       x[2], cd, cg, cb);
  build(oe2, 1.249, 118.4, x[2]+180.0, cd, cg, cb);
}

void Glu::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg,"CG",natom++,name());
  printPDBLine(cd,"CD",natom++,name());
  printPDBLine(oe1,"OE1",natom++,name());
  printPDBLine(oe2,"OE2",natom++,name());
}

void Glu::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg;
  v[natom++]=cd;
  v[natom++]=oe1;
  v[natom++]=oe2;
}

// *** His ******

int His::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG"))
      cg=v;
    else if (!strcmp(name,"ND1"))
      nd1=v;
    else if (!strcmp(name,"CD2"))
      cd2=v;
    else if (!strcmp(name,"CE1"))
      ce1=v;
    else if (!strcmp(name,"NE2"))
      ne2=v;
    else 
      return -1;
  }
  return 0;
}

Vector His::cofm() {
  Vector v;
  v=ca+cb+cg+nd1+cd2+ce1+ne2;
  v/=7.0;
  return v;
}

void His::makeSide() {
  build(cb,  1.530, 110.5,   -120.0, ca,  n,  c);
  build(cg,  1.497, 113.8,       x[0], cb, ca,  n);
  build(nd1, 1.378, 122.7,       x[1], cg, cb, ca);
  build(ce1, 1.345, 109.0,    180.0,nd1, cg, cb);
  build(ne2, 1.319, 111.7,      0.0,ce1, cg, cb);
  build(cd2, 1.319, 106.9,      0.0,ne2,ce1,nd1);
}

void His::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg,"CG",natom++,name());
  printPDBLine(nd1,"ND1",natom++,name());
  printPDBLine(ce1,"CE1",natom++,name());
  printPDBLine(ne2,"NE2",natom++,name());
  printPDBLine(cd2,"CD2",natom++,name());
}

void His::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg;
  v[natom++]=nd1;
  v[natom++]=ce1;
  v[natom++]=ne2;
  v[natom++]=cd2;
}


// *** Ile ******

int Ile::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG1"))
      cg1=v;
    else if (!strcmp(name,"CG2"))
      cg2=v;
    else if (!strcmp(name,"CD") || !strcmp(name,"CD1"))
      cd=v;
    else 
      return -1;
  }
  return 0;
}

Vector Ile::cofm() {
  Vector v;
  v=ca+cb+cg1+cg2+cd;
  v/=5.0;
  return v;
}

void Ile::makeSide() {
  build(cb,  1.540, 111.5,   -120.0, ca,  n,  c);
  build(cg1, 1.530, 110.4,       x[0], cb, ca,  n);
  build(cd,  1.513, 113.8,       x[1],cg1, cb, ca);
  build(cg2, 1.521, 110.5,   x[0]-120, cb, ca,  n);
}

void Ile::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg1,"CG1",natom++,name());
  printPDBLine(cg2,"CG2",natom++,name());
  printPDBLine(cd,"CD",natom++,name());
}

void Ile::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg1;
  v[natom++]=cg2;
  v[natom++]=cd;
}

// *** Leu ******

int Leu::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG"))
      cg=v;
    else if (!strcmp(name,"CD1"))
      cd1=v;
    else if (!strcmp(name,"CD2"))
      cd2=v;
    else 
      return -1;
  }
  return 0;
}

Vector Leu::cofm() {
  Vector v;
  v=ca+cb+cg+cd1+cd2;
  v/=5.0;
  return v;
}

void Leu::makeSide() {
  build(cb,  1.530, 110.5,   -120.0, ca,  n,  c);
  build(cg,  1.530, 116.3,       x[0], cb, ca,  n);
  build(cd1, 1.521, 110.7,       x[1], cg, cb, ca);
  build(cd2, 1.521, 110.7,   x[1]+120, cg, cb, ca);
}

void Leu::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg,"CG",natom++,name());
  printPDBLine(cd1,"CD1",natom++,name());
  printPDBLine(cd2,"CD2",natom++,name());
}

void Leu::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg;
  v[natom++]=cd1;
  v[natom++]=cd2;
}

// *** Lys ******

int Lys::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG"))
      cg=v;
    else if (!strcmp(name,"CD"))
      cd=v;
    else if (!strcmp(name,"CE"))
      ce=v;
    else if (!strcmp(name,"NZ"))
      nz=v;
    else 
      return -1;
  }
  return 0;
}

Vector Lys::cofm() {
  Vector v;
  v=ca+cb+cg+cd+ce+nz;
  v/=6.0;
  return v;
}

void Lys::makeSide() {
  build(cb,  1.530, 110.5, -120.0, ca,  n,  c);
  build(cg,  1.520, 114.1,     x[0], cb, ca,  n);
  build(cd,  1.520, 111.3,     x[1], cg, cb, ca);
  build(ce,  1.520, 111.3,     x[2], cd, cg, cb);
  build(nz,  1.489, 111.9,     x[3], ce, cd, cg);
}

void Lys::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg,"CG",natom++,name());
  printPDBLine(cd,"CD",natom++,name());
  printPDBLine(ce,"CE",natom++,name());
  printPDBLine(nz,"NZ",natom++,name());
}

void Lys::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg;
  v[natom++]=cd;
  v[natom++]=ce;
  v[natom++]=nz;
}


// *** Met ******

int Met::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG"))
      cg=v;
    else if (!strcmp(name,"SD"))
      sd=v;
    else if (!strcmp(name,"CE"))
      ce=v;
    else 
      return -1;
  }
  return 0;
}

Vector Met::cofm() {
  Vector v;
  v=ca+cb+cg+sd+ce;
  v/=5.0;
  return v;
}

void Met::makeSide() {
  build(cb,  1.530, 110.5, -120.0, ca,  n,  c);
  build(cg,  1.520, 114.1,     x[0], cb, ca,  n);
  build(sd,  1.803, 112.7,     x[1], cg, cb, ca);
  build(ce,  1.791, 100.9,     x[2], sd, cg, cb);
}

void Met::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg,"CG",natom++,name());
  printPDBLine(sd,"SD",natom++,name());
  printPDBLine(ce,"CE",natom++,name());
}

void Met::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg;
  v[natom++]=sd;
  v[natom++]=ce;
}


// *** Phe ******

int Phe::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG"))
      cg=v;
    else if (!strcmp(name,"CD1"))
      cd1=v;
    else if (!strcmp(name,"CD2"))
      cd2=v;
    else if (!strcmp(name,"CE1"))
      ce1=v;
    else if (!strcmp(name,"CE2"))
      ce2=v;
    else if (!strcmp(name,"CZ"))
      cz=v;
    else 
      return -1;
  }
  return 0;
}

Vector Phe::cofm() {
  Vector v;
  v=ca+cb+cg+cd1+cd2+ce1+ce2+cz;
  v/=8.0;
  return v;
}

void Phe::makeSide() {
  build(cb,  1.530, 110.5,   -120.0, ca,  n,  c);
  build(cg,  1.502, 113.8,       x[0], cb, ca,  n);
  build(cd1, 1.384, 120.7,       x[1], cg, cb, ca);
  build(cd2, 1.384, 120.7, x[1]-180.0, cg, cb, ca);
  build(ce1, 1.382, 120.7,    180.0,cd1, cg, cb);
  build(ce2, 1.382, 120.7,    180.0,cd2, cg, cb);
  build(cz,  1.382, 120.0,      0.0,ce1,cd1, cg);  
}

void Phe::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg,"CG",natom++,name());
  printPDBLine(cd1,"CD1",natom++,name());
  printPDBLine(cd2,"CD2",natom++,name());
  printPDBLine(ce1,"CE1",natom++,name());
  printPDBLine(ce2,"CE2",natom++,name());
  printPDBLine(cz,"CZ",natom++,name());
}

void Phe::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg;
  v[natom++]=cd1;
  v[natom++]=cd2;
  v[natom++]=ce1;
  v[natom++]=ce2;
  v[natom++]=cz;
}

// *** Pro ******

int Pro::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG"))
      cg=v;
    else if (!strcmp(name,"CD"))
      cd=v;
    else 
      return -1;
  }
  return 0;
}

Vector Pro::cofm() {
  Vector v;
  v=ca+cb+cg+cd;
  v/=4.0;
  return v;
}

void Pro::makeSide() {
  build(cb,  1.530, 104.0,   -120.0, ca,  n,  c);
  build(cg,  1.492, 104.5,       x[0], cb, ca,  n);
  build(cd,  1.503, 106.1,       x[1], cg, cb, ca);
}

void Pro::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg,"CG",natom++,name());
  printPDBLine(cd,"CD",natom++,name());
}

void Pro::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg;
  v[natom++]=cd;
}

// *** Ser ******

int Ser::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"OG"))
      og=v;
    else 
      return -1;
  }
  return 0;
}

Vector Ser::cofm() {
  Vector v;
  v=ca+cb+og;
  v/=3.0;
  return v;
}

void Ser::makeSide() {
  build(cb,  1.530, 110.5,   -120.0, ca,  n,  c);
  build(og,  1.417, 111.1,       x[0], cb, ca,  n);
}

void Ser::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(og,"OG",natom++,name());
}

void Ser::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=og;
}

// *** Thr ******

int Thr::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"OG1"))
      og1=v;
    else if (!strcmp(name,"CG2"))
      cg2=v;
    else 
      return -1;
  }
  return 0;
}

Vector Thr::cofm() {
  Vector v;
  v=ca+cb+og1+cg2;
  v/=4.0;
  return v;
}

void Thr::makeSide() {
  build(cb,  1.540, 111.5,   -120.0, ca,  n,  c);
  build(og1, 1.433, 109.6,       x[0], cb, ca,  n);
  build(cg2, 1.521, 110.5, x[0]-120.0, cb, ca,  n);
}

void Thr::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(og1,"OG1",natom++,name());
  printPDBLine(cg2,"CG2",natom++,name());
}

void Thr::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=og1;
  v[natom++]=cg2;
}


// *** Trp ******

int Trp::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG"))
      cg=v;
    else if (!strcmp(name,"CD1"))
      cd1=v;
    else if (!strcmp(name,"CD2"))
      cd2=v;
    else if (!strcmp(name,"NE1"))
      ne1=v;
    else if (!strcmp(name,"CE2"))
      ce2=v;
    else if (!strcmp(name,"CE3"))
      ce3=v;
    else if (!strcmp(name,"CZ2"))
      cz2=v;
    else if (!strcmp(name,"CZ3"))
      cz3=v;
    else if (!strcmp(name,"CH2"))
      ch2=v;
    else 
      return -1;
  }
  return 0;
}

Vector Trp::cofm() {
  Vector v;
  v=ca+cb+cg+cd1+cd2+ne1+ce2+ce3+cz2+cz3+ch2;
  v/=11.0;
  return v;
}

void Trp::makeSide() {
  build(cb,  1.530, 110.5,   -120.0, ca,  n,  c);
  build(cg,  1.498, 113.6,       x[0], cb, ca,  n);
  build(cd1, 1.365, 126.9,       x[1], cg, cb, ca);
  build(cd2, 1.433, 126.8, x[1]-180.0, cg, cb, ca);
  build(ne1, 1.374, 110.2,    180.0,cd1, cg, cb);
  build(ce2, 1.409, 107.2,    180.0,cd2, cg, cb);
  build(ce3, 1.398, 118.8,    180.0,cd2,ce2,ne1);
  build(cz2, 1.394, 122.4,      0.0,ce2,cd2,ce3);
  build(cz3, 1.382, 118.6,      0.0,ce3,cd2,ce2);
  build(ch2, 1.368, 117.5,      0.0,cz2,ce2,cd2);
}

void Trp::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg,"CG",natom++,name());
  printPDBLine(cd1,"CD1",natom++,name());
  printPDBLine(cd2,"CD2",natom++,name());
  printPDBLine(ne1,"NE1",natom++,name());
  printPDBLine(ce2,"CE2",natom++,name());
  printPDBLine(ce3,"CE3",natom++,name());
  printPDBLine(cz2,"CZ2",natom++,name());
  printPDBLine(cz3,"CZ3",natom++,name());
  printPDBLine(ch2,"CH2",natom++,name());
}

void Trp::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg;
  v[natom++]=cd1;
  v[natom++]=cd2;
  v[natom++]=ne1;
  v[natom++]=ce2;
  v[natom++]=ce3;
  v[natom++]=cz2;
  v[natom++]=cz3;
  v[natom++]=ch2;
}

// *** Tyr ******

int Tyr::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG"))
      cg=v;
    else if (!strcmp(name,"CD1"))
      cd1=v;
    else if (!strcmp(name,"CD2"))
      cd2=v;
    else if (!strcmp(name,"CE1"))
      ce1=v;
    else if (!strcmp(name,"CE2"))
      ce2=v;
    else if (!strcmp(name,"CZ"))
      cz=v;
    else if (!strcmp(name,"OH"))
      oh=v;
    else 
      return -1;
  }
  return 0;
}

Vector Tyr::cofm() {
  Vector v;
  v=ca+cb+cg+cd1+cd2+ce1+ce2+cz+oh;
  v/=9.0;
  return v;
}

void Tyr::makeSide() {
  build(cb,  1.530, 110.5,   -120.0, ca,  n,  c);
  build(cg,  1.512, 113.9,       x[0], cb, ca,  n);
  build(cd1, 1.389, 120.8,       x[1], cg, cb, ca);
  build(cd2, 1.389, 120.8, x[1]-180.0, cg, cb, ca);
  build(ce1, 1.382, 121.2,    180.0,cd1, cg, cb);
  build(ce2, 1.382, 121.2,    180.0,cd2, cg, cb);
  build(cz,  1.378, 119.6,      0.0,ce1,cd1, cg);  
  build(oh,  1.376, 119.9,    180.0, cz,ce2,cd2);  
}
void Tyr::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg,"CG",natom++,name());
  printPDBLine(cd1,"CD1",natom++,name());
  printPDBLine(cd2,"CD2",natom++,name());
  printPDBLine(ce1,"CE1",natom++,name());
  printPDBLine(ce2,"CE2",natom++,name());
  printPDBLine(cz,"CZ",natom++,name());
  printPDBLine(oh,"OH",natom++,name());
}

void Tyr::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg;
  v[natom++]=cd1;
  v[natom++]=cd2;
  v[natom++]=ce1;
  v[natom++]=ce2;
  v[natom++]=cz;
  v[natom++]=oh;
}

// *** Val ******

int Val::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else if (!strcmp(name,"CG1"))
      cg1=v;
    else if (!strcmp(name,"CG2"))
      cg2=v;
    else 
      return -1;
  }
  return 0;
}

Vector Val::cofm() {
  Vector v;
  v=ca+cb+cg1+cg2;
  v/=4.0;
  return v;
}

void Val::makeSide() {
  build(cb,  1.540, 111.5,   -120.0, ca,  n,  c);
  build(cg1, 1.521, 110.5,       x[0], cb, ca,  n);
  build(cg2, 1.521, 110.5, x[0]+120.0, cb, ca,  n);
}

void Val::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
  printPDBLine(cg1,"CG1",natom++,name());
  printPDBLine(cg2,"CG2",natom++,name());
}

void Val::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
  v[natom++]=cg1;
  v[natom++]=cg2;
}

// *** Ala ******

int Ala::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    if (!strcmp(name,"CB"))
      cb=v;
    else 
      return -1;
  }
  return 0;
}

Vector Ala::cofm() {
  Vector v;
  v=ca+cb;
  v/=2.0;
  return v;
}

void Ala::makeSide() {
  build(cb,  1.521, 110.4,   -120.0, ca,  n,  c);  
}

void Ala::print( int& natom) {
  printBack(natom);
  printPDBLine(cb,"CB",natom++,name());
}

void Ala::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
  v[natom++]=cb;
}

// *** Gly ******

int Gly::set(char *name, Vector v) {
  if (setBack(name,v)<0) {
    return -1;
  }
  return 0;
}

Vector Gly::cofm() {
  Vector v;
  v=ca;
  v/=1.0;
  return v;
}

void Gly::makeSide() {
}

void Gly::print(int& natom) {
  printBack(natom);
}

void Gly::addVector(Vector *v, int& natom) {
  addVectorBack(v,natom);
}

// ****************************************************************

AA* getAA(char *resname) {
  AA *aa;

  if (!strcmp(resname,"ARG")) 
    aa=new Arg;
  else if (!strcmp(resname,"ASN")) 
    aa=new Asn;
  else if (!strcmp(resname,"ASP")) 
    aa=new Asp;
  else if (!strcmp(resname,"CYS")) 
    aa=new Cys;
  else if (!strcmp(resname,"GLN"))
    aa=new Gln;
  else if (!strcmp(resname,"GLU"))
    aa=new Glu;
  else if (!strcmp(resname,"HIS") || !strcmp(resname,"HSE") || !strcmp(resname,"HSD"))
    aa=new His;
  else if (!strcmp(resname,"ILE"))
    aa=new Ile;
  else if (!strcmp(resname,"LEU"))
    aa=new Leu;
  else if (!strcmp(resname,"LYS"))
    aa=new Lys;
  else if (!strcmp(resname,"MET"))
    aa=new Met;
  else if (!strcmp(resname,"PHE"))
    aa=new Phe;
  else if (!strcmp(resname,"PRO"))
    aa=new Pro;
  else if (!strcmp(resname,"SER"))
    aa=new Ser;
  else if (!strcmp(resname,"THR"))
    aa=new Thr;
  else if (!strcmp(resname,"TRP"))
    aa=new Trp;
  else if (!strcmp(resname,"TYR"))
    aa=new Tyr;
  else if (!strcmp(resname,"VAL"))
    aa=new Val;
  else if (!strcmp(resname,"ALA"))
    aa=new Ala;
  else if (!strcmp(resname,"GLY"))
    aa=new Gly;
  else {
    fprintf(stderr,"unknown residue %s\n",resname);
    exit(1);
  }
 return aa;
}

int getAminoAcidType(char *resname) {
  static char scres[20][4]={"ARG","ASN","ASP","CYS","GLN",
    			    "GLU","HIS","ILE","LEU","LYS",
			    "MET","PHE","PRO","SER","THR",
			    "TRP","TYR","VAL","ALA","GLY"};
  
  char rname[10];
  if (!strcmp(resname,"HSD") || !strcmp(resname,"HSE") || !strcmp(resname,"HSP")) 
     strcpy(rname,"HIS");
  else 
     strcpy(rname,resname);
 
  for (int i=0; i<20; i++) {
    if (!strcmp(scres[i],rname)) 
      return i;
  } 
 
  fprintf(stderr,"cannot find residue %s\n",resname);
  exit(1);
  return -1;
}  

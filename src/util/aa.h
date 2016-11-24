int getAminoAcidType(char *resname);

class AA {
 protected:
  void nmcrosprod(double  x1, double y1, double z1, 
		    double  x2, double y2, double z2, 
		    double& x3, double& y3, double& z3);

  void viewat(double **M, double **invM, 
		double *P1, double *P2, double *P3);


  void mtxmult(double **invmat, Vector& w);

  char nm[5];

  void printPDBLine(Vector& v, char *aname, int natom, char *resname);

 public:
  AA(char *n, int ndih, int nhea) : 
    nx(ndih), ca(0.0,0.0,0.0), c(0.0,0.0,0.0), n(0.0,0.0,0.0),
    o(0.0,0.0,0.0), sicho(0.0,0.0,0.0), modelca(0), nheavy(nhea) {
    strcpy(nm,n);
    numType=getAminoAcidType(nm);
    pdbnum=-1;
    pdbchain=' ';
    *pdbsegid=0;
  }
  
  Vector ca;
  Vector c;
  Vector n;
  Vector o;

  Vector sicho;

  int numType;
  int pdbnum;
  int modelca;

  char pdbchain;
  char pdbsegid[5];

  int nheavy;

  double x[4];
  int nx;

  void setChi(double chi1, double chi2, double chi3, double chi4);

  int setBack(char *name, Vector v);
  virtual int set(char *name, Vector v)=0;
  virtual Vector cofm()=0;
  virtual void makeSide()=0;
  
  char *name() {return nm;}

  void printBack(int& natom);
  virtual void print(int& natom)=0;

  void addVectorBack(Vector *v, int& natom);
  virtual void addVector(Vector *v, int& natom)=0;

  int heavyAtoms() const { return nheavy; }

  void build(Vector& w, double bond, double angle, double dih, 
	     Vector& v1, Vector &v2, Vector &v3);
};

class Arg : public AA {
 public:
  Arg() : AA("ARG",4,7) {}

  Vector cb;
  Vector cg;
  Vector cd;
  Vector ne;
  Vector cz;
  Vector nh1;
  Vector nh2;
  
  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);

  void addVector(Vector *v, int& natom);
};

class Asn : public AA {
 public:
  Asn() : AA("ASN",2,4) {}

  Vector cb;
  Vector cg;
  Vector od1;
  Vector nd2;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);

  void addVector(Vector *v, int& natom);
};

class Asp : public AA {
 public:
  Asp() : AA("ASP",2,4) {}

  Vector cb;
  Vector cg;
  Vector od1;
  Vector od2;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);

  void addVector(Vector *v, int& natom);
};

class Cys : public AA {
 public:
  Cys() : AA("CYS",1,2) {}

  Vector cb;
  Vector sg;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);

  void addVector(Vector *v, int& natom);
};

class Gln : public AA {
 public:
  Gln() : AA("GLN",3,5) {}

  Vector cb;
  Vector cg;
  Vector cd;
  Vector oe1;
  Vector ne2;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);

  void addVector(Vector *v, int& natom);
};

class Glu : public AA {
 public:
  Glu() : AA("GLU",3,5) {}

  Vector cb;
  Vector cg;
  Vector cd;
  Vector oe1;
  Vector oe2;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);

  void addVector(Vector *v, int& natom);
};

class His : public AA {
 public:
  His() : AA("HIS",2,6) {}

  Vector cb;
  Vector cg;
  Vector nd1;
  Vector cd2;
  Vector ce1;
  Vector ne2;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);

  void addVector(Vector *v, int& natom);
};

class Ile : public AA {
 public:
  Ile() : AA("ILE",2,4) {}

  Vector cb;
  Vector cg1;
  Vector cg2;
  Vector cd;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);

  void addVector(Vector *v, int& natom);
};

class Leu : public AA {
 public:
  Leu() : AA("LEU",2,4) {}

  Vector cb;
  Vector cg;
  Vector cd1;
  Vector cd2;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);
  void addVector(Vector *v, int& natom);
};

class Lys : public AA {
 public:
  Lys() : AA("LYS",4,5) {}

  Vector cb;
  Vector cg;
  Vector cd;
  Vector ce;
  Vector nz;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);
  void addVector(Vector *v, int& natom);
};

class Met : public AA {
 public:
  Met() : AA("MET",3,4) {}

  Vector cb;
  Vector cg;
  Vector sd;
  Vector ce;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);
  void addVector(Vector *v, int& natom);
};

class Phe : public AA {
 public:
  Phe() : AA("PHE",2,7) {}

  Vector cb;
  Vector cg;
  Vector cd1;
  Vector cd2;
  Vector ce1;
  Vector ce2;
  Vector cz;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);
  void addVector(Vector *v, int& natom);
};

class Pro : public AA {
 public:
  Pro() : AA("PRO",2,3) {}

  Vector cb;
  Vector cg;
  Vector cd;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);
  void addVector(Vector *v, int& natom);
};

class Ser : public AA {
 public:
  Ser() : AA("SER",1,2) {}

  Vector cb;
  Vector og;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);
  void addVector(Vector *v, int& natom);
};

class Thr : public AA {
 public:
  Thr() : AA("THR",1,3) {}

  Vector cb;
  Vector og1;
  Vector cg2;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);
  void addVector(Vector *v, int& natom);
};

class Trp : public AA {
 public:
  Trp() : AA("TRP",2,10) {}

  Vector cb;
  Vector cg;
  Vector cd1;
  Vector cd2;
  Vector ne1;
  Vector ce2;
  Vector ce3;
  Vector cz2;
  Vector cz3;
  Vector ch2;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);
  void addVector(Vector *v, int& natom);
};

class Tyr : public AA {
 public:
  Tyr() : AA("TYR",2,8) {}

  Vector cb;
  Vector cg;
  Vector cd1;
  Vector cd2;
  Vector ce1;
  Vector ce2;
  Vector cz;
  Vector oh;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);
  void addVector(Vector *v, int& natom);
};

class Val : public AA {
 public:
  Val() : AA("VAL",1,3) {}

  Vector cb;
  Vector cg1;
  Vector cg2;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);
  void addVector(Vector *v, int& natom);
};

class Ala : public AA {
 public:
  Ala() : AA("ALA",0,1) {}

  Vector cb;

  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);
  void addVector(Vector *v, int& natom);
};

class Gly : public AA {
 public:
  Gly() : AA("GLY",0,0) {}
  int set(char *name, Vector v);
  Vector cofm();
  void makeSide();
  void print(int& natom);
  void addVector(Vector *v, int& natom);
};

AA* getAA(char *resname);


  

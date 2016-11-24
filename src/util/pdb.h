// Michael Feig, 1998

#ifndef PDB_H
#define PDB_H

#include <stdio.h>
#include <string.h>
#include "vector.h"

enum SelEnum { CA, CB, CAB, HEAVY, ALL };

// ### PDBEntry ######

class PDBEntry {
 private:
  int  atminx;
  char atmname[6];
  char resname[6];
  int  resnum;
  char alt;
  char chn;
  Vector coor;
  double aux1,aux2;
  char segid[6];
  int selFlag;

  char *substr(char *l, int offset, int len);

public:
  PDBEntry();
  PDBEntry(int ainx, char *aname, char *rname, int rnum, 
	   Vector v, char *seg=0);
  PDBEntry& operator=(const PDBEntry& from);

  int&    atomIndex()       { return atminx; }
  char*   atomName()        { return atmname; }
  char*   residueName()     { return resname; }
  int&    residueNumber()   { return resnum; }
  char&   alternate()       { return alt; }
  char&   chain()           { return chn; }
  Vector& coordinates()     { return coor; }
  double& auxiliary1()      { return aux1; }
  double& auxiliary2()      { return aux2; }
  char*   segmentID()       { return segid; }

  int&    selected()        { return selFlag; }

  char type();

  int read(FILE *fptr, SelEnum selMode=ALL);
  void write(FILE *fptr);
};

#endif /* PDB_H */












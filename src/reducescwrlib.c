// usage:   reducescwrlib fullLibFile reducedLibFile

// 2000, Michael Feig (meikel@scripps.edu)
// The Scripps Research Institute

// *** C++ code ***

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

class LibRecord {
 public:
  char res[4];
  int phi;
  int psi;
  int N;
  int r1;
  int r2;
  int r3;
  int r4;
  float p;
  float chi1;
  float chi2;
  float chi3;
  float chi4;
};

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

static char res[20][4]={"ARG","ASN","ASP","CYS","GLN",
			"GLU","HIS","ILE","LEU","LYS",
			"MET","PHE","PRO","SER","THR",
			"TRP","TYR","VAL","ALA","GLY"};


int getAminoAcidType(char *resname) {
  for (int i=0; i<20; i++) {
    if (!strcmp(res[i],resname)) 
      return i;
  } 
  
  fprintf(stderr,"cannot find residue %s\n",resname);
  exit(1);

  return -1;
}  


int getPhiPsiIndex(int phi, int psi) {
  int iphi=(phi+180)/10;
  int ipsi=(psi+180)/10;
  
  if (iphi>=36) return -1;
  if (ipsi>=36) return -1;

  return iphi*36+ipsi;
}

int main(int argc, char **argv) {
  if (argc<2) {
    fprintf(stderr,"usage: reducescwrlib fullLIB redLIB\n");
    exit(1);
  }

  LibRecord librecord;
  ReducedLibRecord red;

  int cntaa[20];
  int i;
  for (i=0; i<20; i++) 
    cntaa[i]=0;

  FILE *inptr=fopen(argv[1],"rb");
  FILE *outptr=fopen(argv[2],"wb");
  
  while (!feof(inptr)) {
    if (fread(&librecord,sizeof(LibRecord),1,inptr)) {
      red.phipsi=getPhiPsiIndex(librecord.phi,librecord.psi);
      if (red.phipsi>=0) {
	red.aatype=getAminoAcidType(librecord.res);
	red.p=librecord.p;
	red.chi1=librecord.chi1;
	red.chi2=librecord.chi2;
	red.chi3=librecord.chi3;
	red.chi4=librecord.chi4;
	cntaa[red.aatype]++;
	fwrite(&red,sizeof(ReducedLibRecord),1,outptr);
      }
    }
  }
  fclose(inptr);
  fclose(outptr);
  
  //  for (i=0; i<20; i++) {
  //    printf ("%s: %d\n",res[i],cntaa[i]/36/36);
  //  }

  return 0;
}
